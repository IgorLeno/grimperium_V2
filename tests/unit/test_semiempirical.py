"""
Unit tests for SemiempiricalHandler.

Tests cover:
    - Initialization and configuration
    - Method validation
    - Mock calculation interface

Note: Actual PM7 calculations require MOPAC installation
and are tested separately in integration tests.

"""

import pytest

from grimperium.data.semiempirical import SemiempiricalHandler


class TestSemiempiricalHandler:
    """Tests for SemiempiricalHandler class."""

    def test_init_default(self):
        """Test default initialization."""
        handler = SemiempiricalHandler()

        assert handler.method == "PM7"
        assert handler.use_crest is True
        assert handler.mopac_path is None

    def test_init_pm6(self):
        """Test initialization with PM6 method."""
        handler = SemiempiricalHandler(method="PM6")
        assert handler.method == "PM6"

    def test_init_invalid_method(self):
        """Test that invalid method raises ValueError."""
        with pytest.raises(ValueError, match="not supported"):
            SemiempiricalHandler(method="INVALID")

    def test_supported_methods(self):
        """Test that all supported methods are defined."""
        expected_methods = ["PM7", "PM6", "PM6-D3H+", "AM1", "RM1"]

        for method in expected_methods:
            assert method in SemiempiricalHandler.SUPPORTED_METHODS

    def test_repr(self):
        """Test string representation."""
        handler = SemiempiricalHandler(method="PM7", use_crest=True)
        repr_str = repr(handler)

        assert "SemiempiricalHandler" in repr_str
        assert "PM7" in repr_str
        assert "crest=True" in repr_str

    def test_init_no_crest(self):
        """Test initialization without CREST."""
        handler = SemiempiricalHandler(use_crest=False)
        assert handler.use_crest is False

    def test_work_dir_default(self):
        """Test default work directory."""
        handler = SemiempiricalHandler()
        assert str(handler.work_dir) == "mopac_work"


class TestSemiempiricalHandlerMethods:
    """Tests for SemiempiricalHandler calculation methods."""

    @pytest.mark.skip(reason="Calculate single not implemented yet")
    def test_calculate_single(self, sample_smiles):
        """Test single molecule calculation stub."""
        handler = SemiempiricalHandler()
        # This will raise NotImplementedError until implemented
        with pytest.raises(NotImplementedError):
            handler.calculate_single(sample_smiles[0])

    @pytest.mark.skip(reason="Calculate batch not implemented yet")
    def test_calculate_batch(self, sample_smiles):
        """Test batch calculation stub."""
        handler = SemiempiricalHandler()
        with pytest.raises(NotImplementedError):
            handler.calculate_batch(sample_smiles)

    @pytest.mark.skip(reason="Calculate from xyz not implemented yet")
    def test_calculate_from_xyz(self):
        """Test calculation from XYZ coordinates."""
        handler = SemiempiricalHandler()
        xyz = """3
        water
        O  0.0  0.0  0.0
        H  0.96 0.0  0.0
        H -0.24 0.93 0.0
        """
        with pytest.raises(NotImplementedError):
            handler.calculate_from_xyz(xyz)


class TestSemiempiricalHandlerPM7Choice:
    """Tests documenting why PM7 was chosen."""

    def test_pm7_is_default(self):
        """PM7 should be the default method."""
        handler = SemiempiricalHandler()
        assert handler.method == "PM7"

    def test_pm7_recommended_for_enthalpy(self):
        """
        Document that PM7 is recommended for enthalpy calculations.

        PM7 (Stewart, 2013) is the best general-purpose semiempirical
        for enthalpy of formation because:
        - Specifically parametrized for thermochemistry
        - Excellent equilibrium geometries
        - Good radical description
        - Wide element coverage
        """
        # This test documents the design decision
        assert "PM7" in SemiempiricalHandler.SUPPORTED_METHODS
        # PM7 should be first (default choice)
        assert SemiempiricalHandler.SUPPORTED_METHODS[0] == "PM7"
