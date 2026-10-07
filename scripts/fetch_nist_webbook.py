#!/usr/bin/env python3
"""Fetch NIST Chemistry WebBook pages for the experimental ΔHf validation set.

For each molecular formula of the CBS reference (neutral singlets) and of the
current experimental table, one formula search restricted to species with
gas or condensed-phase thermochemistry (``cTG`` and ``cTC`` together are a
union) is fetched, then the compound page of every listed species with the
gas, condensed and phase-change sections (``Mask=7`` = 1 | 2 | 4).

Why the rules below:

* NIST values are used for validation only: the site states "All rights
  reserved" and its robots.txt carries ``Content-Signal: ai-train=no``;
* robots.txt asks for ``Crawl-delay: 5``: one request every
  ``REQUEST_INTERVAL_S`` seconds, single thread, honest User-Agent without
  personal data;
* any 403/429 or bot-challenge page stops the run at once, so nothing is
  retried against an explicit refusal; other failures are retried once and
  three in a row stop the run;
* pages are cached gzip-compressed under ``--cache-dir`` and written
  atomically, so an interrupted run resumes by skipping cached pages;
  "not found" search pages are cached too, because they are answers.

The cache (``data/raw/``) stays out of git. Parsing of compound pages lives in
``build_experimental_hf.py``; this script only extracts compound ids.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import re
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from collections import Counter
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path

BASE_URL = "https://webbook.nist.gov/cgi/cbook.cgi"
USER_AGENT = "grimperium-research/0.1 (academic thermochemistry validation)"
REQUEST_INTERVAL_S = 5.0
RETRY_WAIT_S = 60.0
MAX_CONSECUTIVE_FAILURES = 3
TIMEOUT_S = 60.0
#: Gas | condensed | phase change sections of a compound page.
COMPOUND_MASK = 7
STOP_STATUSES = frozenset({403, 429})
#: Interstitial challenge pages only: every normal page already loads the
#: ``/cdn-cgi/challenge-platform`` detection script, so that is no marker.
CHALLENGE_MARKERS = ("Just a moment...", "cf_chl_opt", "Attention Required!", "captcha")
WEBBOOK_MARKER = "NIST Chemistry WebBook"

RESULT_COUNT = re.compile(r"(\d+)\s+matching species were found")
#: Result item: compound id and the formula shown after the name.
RESULT_ITEM = re.compile(
    r'<li><a href="/cgi/cbook\.cgi\?ID=([A-Za-z0-9]+)&amp;[^"]*">[^<]*</a>'
    r"\s*\(((?:[^()<]|<sub>|</sub>)*)\)"
)
COMPOUND_LINK = re.compile(r"/cgi/cbook\.cgi\?ID=([A-Za-z0-9]+)&amp;Units=SI&amp;Mask=")
NOT_FOUND = "Chemical Formula Not Found"
FORMULA_SAFE = re.compile(r"^[A-Z][A-Za-z0-9]*$")


class StopFetch(RuntimeError):
    """The server refused or challenged us; stop instead of retrying."""


def collect_formulas(cbs_csv: Path, exp_csv: Path | None) -> list[str]:
    """C-containing formulas of neutral singlet CBS molecules and the exp table."""
    import pandas as pd
    from rdkit import Chem, RDLogger
    from rdkit.Chem.rdMolDescriptors import CalcMolFormula

    RDLogger.DisableLog("rdApp.*")  # type: ignore[attr-defined]
    cbs = pd.read_csv(cbs_csv, usecols=["smiles", "multiplicity", "charge"])
    cbs = cbs[(cbs["multiplicity"] == 1) & (cbs["charge"] == 0)]
    formulas: set[str] = set()
    for smiles in cbs["smiles"].unique():
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            formulas.add(CalcMolFormula(mol))
    if exp_csv is not None:
        exp = pd.read_csv(exp_csv, usecols=["formula"])
        formulas.update(str(f) for f in exp["formula"].dropna())
    bad = sorted(f for f in formulas if not FORMULA_SAFE.match(f))
    if bad:
        raise ValueError(f"unexpected formula strings: {bad[:5]}")
    # The experimental table requires carbon; Hill formulas start with C if any.
    return sorted(f for f in formulas if f.startswith("C"))


def search_url(formula: str) -> str:
    query = {"Formula": formula, "NoIon": "on", "cTG": "on", "cTC": "on", "Units": "SI"}
    return f"{BASE_URL}?{urllib.parse.urlencode(query)}"


def compound_url(compound_id: str) -> str:
    query = {"ID": compound_id, "Units": "SI", "Mask": str(COMPOUND_MASK)}
    return f"{BASE_URL}?{urllib.parse.urlencode(query)}"


def compound_ids(page: str, formula: str | None = None) -> tuple[list[str], str]:
    """Compound ids of a formula search page and the page kind.

    With ``formula``, listed species whose shown formula differs (deuterated
    or tritiated isotopologues such as ``CHD3O`` for ``CH4O``) are skipped:
    the experimental table rejects isotopic labels anyway.

    Kinds: ``list`` (result list), ``not_found``, ``single`` (the search
    matched one species and the server answered with its page directly),
    ``truncated`` (fewer items listed than reported) and ``unknown``.
    """
    if NOT_FOUND in page:
        return [], "not_found"
    count = RESULT_COUNT.search(page)
    if count:
        items = RESULT_ITEM.findall(page)
        kind = "list" if len(items) == int(count.group(1)) else "truncated"
        ids = [
            compound_id
            for compound_id, shown in items
            if formula is None or re.sub(r"</?sub>|\s", "", shown) == formula
        ]
        return list(dict.fromkeys(ids)), kind
    if 'class="inchi-text"' in page or "CAS Registry Number" in page:
        # The page links its own other sections; isotopologue links are rare.
        links = Counter(COMPOUND_LINK.findall(page))
        if links:
            return [links.most_common(1)[0][0]], "single"
    return [], "unknown"


class Fetcher:
    def __init__(self, cache_dir: Path, interval: float = REQUEST_INTERVAL_S) -> None:
        self.cache_dir = cache_dir
        self.interval = interval
        self.log_path = cache_dir / "fetch_log.jsonl"
        self._last_request = 0.0
        self.requests = 0
        self.failures_in_row = 0

    def _log(self, record: dict[str, object]) -> None:
        record["time"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        with self.log_path.open("a", encoding="utf-8") as handle:
            handle.write(json.dumps(record) + "\n")

    def _wait_turn(self) -> None:
        delay = self._last_request + self.interval - time.monotonic()
        if delay > 0:
            time.sleep(delay)
        self._last_request = time.monotonic()

    def _request(self, url: str) -> str:
        if not url.startswith(BASE_URL):
            raise ValueError(f"refusing non-WebBook URL {url}")
        request = urllib.request.Request(  # noqa: S310 - fixed HTTPS base URL
            url, headers={"User-Agent": USER_AGENT}
        )
        for attempt in (1, 2):
            self._wait_turn()
            self.requests += 1
            try:
                with urllib.request.urlopen(  # noqa: S310 - checked above
                    request, timeout=TIMEOUT_S
                ) as response:
                    body: bytes = response.read()
                    status = int(response.status)
            except urllib.error.HTTPError as error:
                self._log({"url": url, "status": error.code, "attempt": attempt})
                if error.code in STOP_STATUSES:
                    raise StopFetch(f"HTTP {error.code} for {url}") from error
            except (urllib.error.URLError, TimeoutError, OSError) as error:
                self._log({"url": url, "error": str(error), "attempt": attempt})
            else:
                page = body.decode("utf-8", errors="replace")
                if any(marker in page for marker in CHALLENGE_MARKERS):
                    self._log({"url": url, "status": status, "challenge": True})
                    raise StopFetch(f"bot challenge served for {url}")
                if WEBBOOK_MARKER not in page:
                    self._log({"url": url, "status": status, "not_webbook": True})
                    raise StopFetch(f"unexpected page for {url}")
                self.failures_in_row = 0
                self._log({"url": url, "status": status, "bytes": len(body)})
                return page
            if attempt == 1:
                time.sleep(RETRY_WAIT_S)
        self.failures_in_row += 1
        if self.failures_in_row >= MAX_CONSECUTIVE_FAILURES:
            raise StopFetch(f"{MAX_CONSECUTIVE_FAILURES} failed pages in a row")
        raise OSError(f"failed twice: {url}")

    def get(self, relative: str, url: str) -> tuple[str, bool]:
        """Return (page, from_cache); fetch and cache atomically if missing."""
        path = self.cache_dir / relative
        if path.exists():
            with gzip.open(path, "rt", encoding="utf-8") as handle:
                return handle.read(), True
        page = self._request(url)
        path.parent.mkdir(parents=True, exist_ok=True)
        partial = path.with_name(path.name + ".part")
        with gzip.open(partial, "wt", encoding="utf-8") as handle:
            handle.write(page)
        os.replace(partial, path)
        return page, False


def run(fetcher: Fetcher, formulas: Iterable[str]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    formula_list = list(formulas)
    for index, formula in enumerate(formula_list, start=1):
        try:
            page, cached = fetcher.get(f"search/{formula}.html.gz", search_url(formula))
        except OSError as error:
            counts["search_failed"] += 1
            print(f"[{index}/{len(formula_list)}] {formula}: {error}", flush=True)
            continue
        ids, kind = compound_ids(page, formula)
        counts[f"search_{kind}"] += 1
        if kind in {"truncated", "unknown"}:
            print(f"warning: {formula} search page is {kind}", flush=True)
        fetched = 0
        for compound_id in ids:
            relative = f"compound/{compound_id}.html.gz"
            try:
                _, compound_cached = fetcher.get(relative, compound_url(compound_id))
            except OSError as error:
                counts["compound_failed"] += 1
                print(f"  {compound_id}: {error}", flush=True)
                continue
            fetched += int(not compound_cached)
        counts["compound_ids"] += len(ids)
        print(
            f"[{index}/{len(formula_list)}] {formula}: {kind}, {len(ids)} ids, "
            f"{fetched} fetched, search {'cached' if cached else 'fetched'}, "
            f"requests {fetcher.requests}",
            flush=True,
        )
    return dict(counts)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cbs-csv", type=Path, required=True)
    parser.add_argument("--exp-csv", type=Path, default=None)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=None, help="first N formulas")
    parser.add_argument("--formulas", nargs="*", default=None, help="explicit list")
    args = parser.parse_args(argv)

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    formulas = args.formulas or collect_formulas(args.cbs_csv, args.exp_csv)
    if not args.formulas:
        (args.cache_dir / "formulas.txt").write_text("\n".join(formulas) + "\n")
    if args.limit is not None:
        formulas = formulas[: args.limit]
    print(f"{len(formulas)} formulas; cache {args.cache_dir}", flush=True)

    fetcher = Fetcher(args.cache_dir)
    try:
        counts = run(fetcher, formulas)
    except StopFetch as stop:
        print(f"STOPPED: {stop}", file=sys.stderr, flush=True)
        return 2
    print(json.dumps({"requests": fetcher.requests, **counts}), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
