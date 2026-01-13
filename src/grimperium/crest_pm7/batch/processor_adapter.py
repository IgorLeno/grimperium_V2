"""Processor adapter for batch processing with fixed timeouts.

This module provides:
- FixedTimeoutPredictor: Returns fixed timeouts instead of dynamic predictions
- FixedTimeoutProcessor: Wraps MoleculeProcessor for batch processing

Used by BatchExecutionManager to process molecules with batch-level timeouts.
"""

import logging
from collections import deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from grimperium.crest_pm7.config import PM7Config, TimeoutConfidence
from grimperium.crest_pm7.molecule_processor import MoleculeProcessor, PM7Result

LOG = logging.getLogger("grimperium.crest_pm7.batch.processor_adapter")


@dataclass
class FixedTimeoutPredictor:
    """A timeout predictor that always returns fixed values.

    Used by FixedTimeoutProcessor to override dynamic timeout prediction
    with batch-level fixed timeouts.

    Implements the same interface as TimeoutPredictor for compatibility.

    Attributes:
        crest_timeout_seconds: Fixed CREST timeout in seconds
        mopac_timeout_seconds: Fixed MOPAC timeout in seconds
        observations: Bounded deque of (nheavy, time) tuples for statistics
        max_observations: Maximum observations to keep (default 100)
    """

    crest_timeout_seconds: float
    mopac_timeout_seconds: float
    max_observations: int = 100
    observations: deque = field(default_factory=lambda: deque(maxlen=100))

    def __post_init__(self):
        """Initialize deque with proper maxlen."""
        if not isinstance(self.observations, deque):
            self.observations = deque(self.observations, maxlen=self.max_observations)
        elif self.observations.maxlen != self.max_observations:
            # Rebuild deque with correct maxlen
            self.observations = deque(self.observations, maxlen=self.max_observations)

    # Properties for TimeoutPredictor interface compatibility
    @property
    def n_samples(self) -> int:
        """Number of observations collected."""
        return len(self.observations)

    @property
    def is_fitted(self) -> bool:
        """Always returns True (no fitting needed for fixed timeout)."""
        return True

    def predict(
        self,
        nheavy: int,  # noqa: ARG002 - required for interface compatibility
        num_conformers: int = 1,  # noqa: ARG002 - required for interface compatibility
    ) -> tuple[float, TimeoutConfidence]:
        """Return fixed timeout regardless of molecule complexity.

        Args:
            nheavy: Number of heavy atoms (ignored)
            num_conformers: Number of conformers (ignored)

        Returns:
            Tuple of (mopac_timeout_seconds, HIGH confidence)
        """
        # For batch processing, we return the total MOPAC timeout
        # The processor will distribute it across conformers
        return self.mopac_timeout_seconds, TimeoutConfidence.HIGH

    def add_observation(self, nheavy: int, execution_time: float) -> None:
        """Record observation for statistics (no model update).

        Automatically evicts old observations if maxlen exceeded.

        Args:
            nheavy: Number of heavy atoms
            execution_time: Actual execution time
        """
        self.observations.append((nheavy, execution_time))
        LOG.debug(f"Recorded observation: nheavy={nheavy}, time={execution_time:.1f}s")

    def clear_observations(self) -> None:
        """Clear observations for new batch."""
        self.observations.clear()
        LOG.debug("Cleared timeout observations")

    def fit(self) -> bool:
        """No-op for fixed timeout predictor.

        Returns:
            Always True
        """
        return True

    def get_stats(self) -> dict:
        """Get predictor statistics.

        Returns:
            Dict with fixed timeout values and observation stats
        """
        stats = {
            "type": "fixed",
            "crest_timeout_seconds": self.crest_timeout_seconds,
            "mopac_timeout_seconds": self.mopac_timeout_seconds,
            "n_observations": self.n_samples,
        }

        if self.observations:
            times = [t for _, t in self.observations]
            stats["observed_time_range"] = (min(times), max(times))
            stats["observed_time_avg"] = sum(times) / len(times)

        return stats


class FixedTimeoutProcessor:
    """Wraps MoleculeProcessor for batch processing with fixed timeouts.

    This adapter:
    1. Creates a MoleculeProcessor with FixedTimeoutPredictor
    2. Provides a simpler interface for batch processing
    3. Tracks actual timeouts used

    Attributes:
        config: PM7 configuration
        crest_timeout_minutes: Fixed CREST timeout
        mopac_timeout_minutes: Fixed MOPAC timeout
        processor: Underlying MoleculeProcessor
    """

    def __init__(
        self,
        config: PM7Config,
        crest_timeout_minutes: float = 30.0,
        mopac_timeout_minutes: float = 60.0,
    ) -> None:
        """Initialize processor with fixed timeouts.

        Args:
            config: PM7 configuration
            crest_timeout_minutes: Fixed CREST timeout in minutes
            mopac_timeout_minutes: Fixed MOPAC timeout in minutes
        """
        self.config = config
        self.crest_timeout_minutes = crest_timeout_minutes
        self.mopac_timeout_minutes = mopac_timeout_minutes

        # Create fixed timeout predictor
        self._timeout_predictor = FixedTimeoutPredictor(
            crest_timeout_seconds=crest_timeout_minutes * 60,
            mopac_timeout_seconds=mopac_timeout_minutes * 60,
        )

        # Create processor with fixed timeout predictor
        # Note: FixedTimeoutPredictor implements TimeoutPredictor protocol
        self.processor = MoleculeProcessor(
            config=config,
            timeout_predictor=self._timeout_predictor,  # type: ignore[arg-type]
        )

        LOG.info(
            f"FixedTimeoutProcessor initialized: "
            f"CREST={crest_timeout_minutes}min, MOPAC={mopac_timeout_minutes}min"
        )

    def process_with_fixed_timeout(
        self,
        mol_id: str,
        smiles: str,
        input_xyz: Optional[Path] = None,
    ) -> PM7Result:
        """Process molecule with fixed timeouts.

        This is the main entry point for batch processing.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            input_xyz: Optional input XYZ file

        Returns:
            PM7Result from processing
        """
        LOG.debug(f"Processing {mol_id} with fixed timeouts")
        return self.processor.process(mol_id, smiles, input_xyz)

    def update_timeouts(
        self,
        crest_timeout_minutes: float,
        mopac_timeout_minutes: float,
    ) -> None:
        """Update timeout values for next batch.

        Args:
            crest_timeout_minutes: New CREST timeout
            mopac_timeout_minutes: New MOPAC timeout
        """
        self.crest_timeout_minutes = crest_timeout_minutes
        self.mopac_timeout_minutes = mopac_timeout_minutes
        self._timeout_predictor.crest_timeout_seconds = crest_timeout_minutes * 60
        self._timeout_predictor.mopac_timeout_seconds = mopac_timeout_minutes * 60

        LOG.debug(
            f"Updated timeouts: CREST={crest_timeout_minutes}min, "
            f"MOPAC={mopac_timeout_minutes}min"
        )

    def get_timeout_stats(self) -> dict:
        """Get timeout statistics from observations.

        Returns:
            Dict with timeout configuration and observation stats
        """
        return self._timeout_predictor.get_stats()
