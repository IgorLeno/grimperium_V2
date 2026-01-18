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

from grimperium import DequeAny, DictStrAny
from grimperium.crest_pm7.config import PM7Config, TimeoutConfidence
from grimperium.crest_pm7.molecule_processor import MoleculeProcessor, PM7Result
from grimperium.crest_pm7.preoptimization import xTBPreOptimizer

LOG = logging.getLogger("grimperium.crest_pm7.batch.processor_adapter")

# Default maximum observations to keep for statistics
DEFAULT_MAX_OBSERVATIONS = 100


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
    max_observations: int = DEFAULT_MAX_OBSERVATIONS
    # Use unbounded deque as placeholder, __post_init__ creates real one
    observations: DequeAny = field(default_factory=deque)

    def __post_init__(self) -> None:
        """Initialize observations deque with correct maxlen.

        This is called after dataclass __init__, so max_observations
        is already set to the provided value (or default).
        """
        self.observations = deque(maxlen=self.max_observations)
        LOG.debug(f"Initialized observations deque with maxlen={self.max_observations}")

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

    def get_stats(self) -> DictStrAny:
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
    4. Optionally runs xTB pre-optimization before CREST

    Attributes:
        config: PM7 configuration
        crest_timeout_minutes: Fixed CREST timeout
        mopac_timeout_minutes: Fixed MOPAC timeout
        processor: Underlying MoleculeProcessor
        preoptimizer: xTB pre-optimizer (enabled via constructor)
    """

    def __init__(
        self,
        config: PM7Config,
        crest_timeout_minutes: float = 30.0,
        mopac_timeout_minutes: float = 60.0,
        enable_xtb_preopt: bool = False,
    ) -> None:
        """Initialize FixedTimeoutProcessor with batch-specific timeout.

        CRITICAL: This adapter mutates the passed config object's crest_timeout
        field to enforce the batch-specific timeout. This is intentional:
        run_crest() uses config.crest_timeout for subprocess timeout values.

        If the config object is shared with other components or threads,
        this mutation could cause unexpected behavior. To prevent this,
        pass a copy of the config:

            from copy import copy
            config_copy = copy(config)
            adapter = ProcessorAdapter(config_copy, timeout_minutes=20)

        Args:
            config: PM7Config instance (will be mutated)
            crest_timeout_minutes: Fixed CREST timeout in minutes
            mopac_timeout_minutes: Fixed MOPAC timeout in minutes
            enable_xtb_preopt: Enable xTB pre-optimization (default: False)

        Side Effects:
            Modifies config.crest_timeout = crest_timeout_minutes * 60
            This ensures run_crest() respects the batch timeout.
        """
        self.config = config
        self.crest_timeout_minutes = crest_timeout_minutes
        self.mopac_timeout_minutes = mopac_timeout_minutes

        # CRITICAL: Update config.crest_timeout to match batch timeout
        # run_crest() uses config.crest_timeout for subprocess timeout
        self.config.crest_timeout = crest_timeout_minutes * 60  # Convert to seconds

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

        # Create xTB pre-optimizer
        self.preoptimizer = xTBPreOptimizer(enabled=enable_xtb_preopt)

        LOG.info(
            f"FixedTimeoutProcessor initialized: "
            f"CREST={crest_timeout_minutes}min, MOPAC={mopac_timeout_minutes}min, "
            f"xTB_preopt={enable_xtb_preopt}"
        )

    def process_with_fixed_timeout(
        self,
        mol_id: str,
        smiles: str,
        input_xyz: Path | None = None,
    ) -> PM7Result:
        """Process molecule with fixed timeouts.

        This is the main entry point for batch processing.
        If xTB pre-optimization is enabled and input_xyz is provided,
        the structure will be pre-optimized before CREST processing.

        Args:
            mol_id: Molecule identifier
            smiles: SMILES string
            input_xyz: Optional input XYZ file

        Returns:
            PM7Result from processing
        """
        LOG.debug(f"Processing {mol_id} with fixed timeouts")

        # xTB pre-optimization before MoleculeProcessor
        if input_xyz is not None and self.preoptimizer.enabled:
            work_dir = self.config.temp_dir / mol_id
            work_dir.mkdir(parents=True, exist_ok=True)

            preopt_result = self.preoptimizer.preoptimize_structure(
                mol_id, input_xyz, work_dir
            )
            if preopt_result.success and preopt_result.output_xyz:
                input_xyz = preopt_result.output_xyz
                LOG.info(f"Using xTB pre-optimized structure: {input_xyz}")
            else:
                LOG.warning(f"xTB pre-opt failed for {mol_id}, using original")

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

        # CRITICAL: Also update config.crest_timeout for run_crest()
        # This ensures the subprocess.run() timeout is correct
        self.config.crest_timeout = crest_timeout_minutes * 60  # Convert to seconds

        LOG.debug(
            f"Updated timeouts: CREST={crest_timeout_minutes}min, "
            f"MOPAC={mopac_timeout_minutes}min"
        )

    def get_timeout_stats(self) -> DictStrAny:
        """Get timeout statistics from observations.

        Returns:
            Dict with timeout configuration and observation stats
        """
        return self._timeout_predictor.get_stats()
