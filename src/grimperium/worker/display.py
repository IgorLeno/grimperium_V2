"""Live molecule progress display for the worker process.

Renders a one-row Rich Live table showing the molecule currently being
processed. Designed for single-molecule-at-a-time workers; not connected
to the CSV-polling ProgressTracker used by the main CLI batch runner.
"""

import time
from dataclasses import dataclass, field

from rich.console import Console
from rich.live import Live
from rich.table import Table


@dataclass
class MoleculeProgressDisplay:
    """Minimal Rich Live display for a single molecule in flight.

    Columns: Mol ID | SMILES (truncated to 30 chars) | Status | Elapsed | Result

    Usage::

        display = MoleculeProgressDisplay()
        display.start()
        display.update("mol_001", "CCO", "CREST")
        # ... processing ...
        display.update("mol_001", "CCO", "DONE", result="OK")
        display.stop()
    """

    console: Console = field(default_factory=Console)
    _live: Live | None = field(default=None, init=False, repr=False)
    _mol_id: str = field(default="", init=False, repr=False)
    _smiles: str = field(default="", init=False, repr=False)
    _status: str = field(default="IDLE", init=False, repr=False)
    _result: str = field(default="", init=False, repr=False)
    _start_time: float = field(default=0.0, init=False, repr=False)

    def _build_table(self) -> Table:
        elapsed = time.monotonic() - self._start_time if self._start_time else 0.0
        elapsed_str = _format_elapsed(elapsed)
        smiles_short = (self._smiles[:27] + "...") if len(self._smiles) > 30 else self._smiles

        table = Table(show_header=True, header_style="bold", box=None, padding=(0, 1))
        table.add_column("Mol ID", min_width=10)
        table.add_column("SMILES", min_width=32)
        table.add_column("Status", min_width=10)
        table.add_column("Elapsed", min_width=8, justify="right")
        table.add_column("Result", min_width=8)

        status_style = _status_style(self._status)
        table.add_row(
            self._mol_id or "—",
            smiles_short or "—",
            f"[{status_style}]{self._status}[/{status_style}]",
            elapsed_str,
            self._result or "—",
        )
        return table

    def start(self) -> None:
        """Start the Live display. Safe to call multiple times."""
        if self._live is None:
            self._live = Live(
                self._build_table(),
                console=self.console,
                refresh_per_second=4,
                transient=True,
            )
            self._live.start()

    def stop(self) -> None:
        """Stop the Live display and clear it from the terminal."""
        if self._live is not None:
            self._live.stop()
            self._live = None

    def update(
        self,
        mol_id: str,
        smiles: str,
        status: str,
        result: str = "",
    ) -> None:
        """Update the displayed state.

        Args:
            mol_id: Molecule identifier.
            smiles: SMILES string (truncated to 30 chars in display).
            status: Current stage label, e.g. "CREST", "MOPAC", "DONE", "FAILED".
            result: Optional outcome label shown in the Result column.
        """
        if self._mol_id != mol_id:
            self._start_time = time.monotonic()

        self._mol_id = mol_id
        self._smiles = smiles
        self._status = status
        self._result = result

        if self._live is not None:
            self._live.update(self._build_table())

    def render_snapshot(self) -> str:
        """Return a plain-text snapshot of the current table for testing."""
        from io import StringIO
        buf = StringIO()
        capture = Console(file=buf, highlight=False, markup=False)
        capture.print(self._build_table())
        return buf.getvalue()


def _format_elapsed(seconds: float) -> str:
    """Format elapsed seconds as HH:MM:SS."""
    total = int(seconds)
    h, remainder = divmod(total, 3600)
    m, s = divmod(remainder, 60)
    if h:
        return f"{h:02d}:{m:02d}:{s:02d}"
    return f"{m:02d}:{s:02d}"


def _status_style(status: str) -> str:
    """Map status label to a Rich style colour."""
    upper = status.upper()
    if upper in ("DONE", "OK", "SUCCESS"):
        return "green"
    if upper in ("FAILED", "ERROR", "SKIP"):
        return "red"
    if upper in ("CREST", "MOPAC", "RUNNING"):
        return "yellow"
    if upper == "IDLE":
        return "dim"
    return "cyan"
