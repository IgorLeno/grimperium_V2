"""Search-page parsing, caching and stop rules of the NIST WebBook fetcher."""

from __future__ import annotations

import gzip
import importlib.util
import io
import sys
import urllib.error
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest


def _load_fetcher() -> ModuleType:
    repo_root = Path(__file__).resolve().parents[2]
    script_path = repo_root / "scripts" / "fetch_nist_webbook.py"
    spec = importlib.util.spec_from_file_location("fetch_nist_webbook", script_path)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


fetcher_module = _load_fetcher()

# Synthetic pages in the WebBook layout (not NIST data).
HEAD = "<html><body><h1>NIST Chemistry WebBook, SRD 69</h1><main>"
LIST_PAGE = (
    HEAD + "<h1>Search Results</h1><p> 3 matching species were found. </p><ol>"
    '<li><a href="/cgi/cbook.cgi?ID=C111&amp;Units=SI&amp;Mask=1">A</a>  '
    "(CH<sub>4</sub>O)<br /></li>"
    '<li><a href="/cgi/cbook.cgi?ID=U222&amp;Units=SI&amp;Mask=2">B</a>  '
    "(CHD<sub>3</sub>O)</li>"
    '<li><a href="/cgi/cbook.cgi?ID=C444&amp;Units=SI&amp;Mask=2">X[0.0<sup>3,6'
    "</sup>]</a>  (CH<sub>4</sub>O)</li>"
    "</ol></main></body></html>"
)
SINGLE_PAGE = (
    HEAD + '<h1 id="Top">A</h1><span class="inchi-text">InChI=1S/CH4/h1H4</span>'
    '<a href="/cgi/cbook.cgi?ID=C333&amp;Units=SI&amp;Mask=2">x</a>'
    '<a href="/cgi/cbook.cgi?ID=C333&amp;Units=SI&amp;Mask=4">y</a>'
    '<a href="/cgi/cbook.cgi?ID=C999&amp;Units=SI&amp;Mask=1">isotopologue</a>'
    "</main></body></html>"
)
NOT_FOUND_PAGE = HEAD + "<h1>Chemical Formula Not Found</h1></main></body></html>"


@pytest.mark.parametrize(
    ("page", "expected"),
    [
        (LIST_PAGE, (["C111", "U222", "C444"], "list")),
        (
            LIST_PAGE.replace(" 3 matching", " 4 matching"),
            (["C111", "U222", "C444"], "truncated"),
        ),
        (SINGLE_PAGE, (["C333"], "single")),
        (
            HEAD + '<span class="inchi-text">InChI=1S/CH4/h1H4</span>'
            '<a href="/cgi/cbook.cgi?Str2File=C555">2d Mol file</a></main>',
            (["C555"], "single"),
        ),
        (NOT_FOUND_PAGE, ([], "not_found")),
        (HEAD + "</main>", ([], "unknown")),
    ],
)
def test_compound_ids(page: str, expected: tuple[list[str], str]) -> None:
    assert fetcher_module.compound_ids(page) == expected


def test_search_url_requests_gas_or_condensed_without_ions() -> None:
    url = fetcher_module.search_url("C2H6O")
    assert "Formula=C2H6O" in url
    assert all(flag in url for flag in ("NoIon=on", "cTG=on", "cTC=on", "Units=SI"))
    assert fetcher_module.compound_url("C64175").endswith("ID=C64175&Units=SI&Mask=7")


class _Response(io.BytesIO):
    status = 200

    def __enter__(self) -> _Response:
        return self

    def __exit__(self, *args: object) -> None:
        self.close()


def _fake_urlopen(pages: list[Any]) -> Any:
    calls: list[str] = []

    def urlopen(request: Any, timeout: float) -> _Response:
        calls.append(request.full_url)
        assert "@" not in str(request.headers)  # no personal data in headers
        item = pages.pop(0)
        if isinstance(item, Exception):
            raise item
        return _Response(item.encode())

    urlopen.calls = calls  # type: ignore[attr-defined]
    return urlopen


def _http_error(code: int) -> urllib.error.HTTPError:
    return urllib.error.HTTPError("u", code, "x", {}, None)  # type: ignore[arg-type]


def test_get_caches_atomically_and_resumes(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    urlopen = _fake_urlopen([LIST_PAGE])
    monkeypatch.setattr(fetcher_module.urllib.request, "urlopen", urlopen)
    fetcher = fetcher_module.Fetcher(tmp_path, interval=0.0)

    page, cached = fetcher.get("search/C2H6O.html.gz", fetcher_module.BASE_URL + "/1")
    assert (page, cached) == (LIST_PAGE, False)
    with gzip.open(tmp_path / "search" / "C2H6O.html.gz", "rt") as handle:
        assert handle.read() == LIST_PAGE
    assert not list(tmp_path.rglob("*.part"))

    again = fetcher_module.Fetcher(tmp_path, interval=0.0)
    assert again.get("search/C2H6O.html.gz", fetcher_module.BASE_URL + "/1") == (
        LIST_PAGE,
        True,
    )
    assert urlopen.calls == [fetcher_module.BASE_URL + "/1"]
    assert (tmp_path / "fetch_log.jsonl").read_text().count('"status": 200') == 1


@pytest.mark.parametrize(
    "response",
    [_http_error(403), _http_error(429), "<title>Just a moment...</title>", "plain"],
)
def test_refusal_or_challenge_stops_without_retry(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, response: Any
) -> None:
    urlopen = _fake_urlopen([response, LIST_PAGE])
    monkeypatch.setattr(fetcher_module.urllib.request, "urlopen", urlopen)
    fetcher = fetcher_module.Fetcher(tmp_path, interval=0.0)
    with pytest.raises(fetcher_module.StopFetch):
        fetcher.get("search/X.html.gz", fetcher_module.BASE_URL + "/1")
    assert len(urlopen.calls) == 1
    assert not (tmp_path / "search" / "X.html.gz").exists()


def test_normal_page_with_detection_script_is_not_a_challenge(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    page = LIST_PAGE + "<script>a.src='/cdn-cgi/challenge-platform/main.js'</script>"
    monkeypatch.setattr(fetcher_module.urllib.request, "urlopen", _fake_urlopen([page]))
    fetcher = fetcher_module.Fetcher(tmp_path, interval=0.0)
    assert fetcher.get("search/X.html.gz", fetcher_module.BASE_URL + "/1") == (
        page,
        False,
    )


def test_transient_errors_retry_once_then_stop_after_three(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(fetcher_module, "RETRY_WAIT_S", 0.0)
    urlopen = _fake_urlopen([_http_error(503), LIST_PAGE] + [_http_error(500)] * 6)
    monkeypatch.setattr(fetcher_module.urllib.request, "urlopen", urlopen)
    fetcher = fetcher_module.Fetcher(tmp_path, interval=0.0)

    assert fetcher.get("search/A.html.gz", fetcher_module.BASE_URL + "/a") == (
        LIST_PAGE,
        False,
    )
    for name in ("B", "C"):
        with pytest.raises(OSError):
            fetcher.get(f"search/{name}.html.gz", f"{fetcher_module.BASE_URL}/{name}")
    with pytest.raises(fetcher_module.StopFetch):
        fetcher.get("search/D.html.gz", fetcher_module.BASE_URL + "/d")
    assert len(urlopen.calls) == 8


def test_isotopologues_are_skipped_by_listed_formula() -> None:
    assert fetcher_module.compound_ids(LIST_PAGE, "CH4O") == (["C111", "C444"], "list")
