"""Tests for configuration schema and microscope config loading."""

import pytest

from merfish_pipeline.config.loader import load_microscope_config
from merfish_pipeline.config.schema import MicroscopeConfig


class TestMicroscopeConfigLoading:
    """Verify all three microscope configs load without errors."""

    @pytest.mark.parametrize("name", ["oni", "nikon", "andor"])
    def test_load_microscope_config(self, name: str):
        cfg = load_microscope_config(name)
        assert isinstance(cfg, MicroscopeConfig)
        assert cfg.file_pattern
        assert cfg.microns_per_pixel > 0
        assert len(cfg.image_dimensions) == 2

    def test_oni_defaults(self):
        cfg = load_microscope_config("oni")
        assert cfg.default_bead_channel == "488nm, Raw"
        assert cfg.position_format == "csv"
        assert cfg.position_file_pattern is not None

    def test_nikon_defaults(self):
        cfg = load_microscope_config("nikon")
        assert cfg.default_bead_channel == "473nm, Raw"
        assert cfg.flip_horizontal is True
        assert cfg.flip_vertical is True

    def test_andor_nullable_fields(self):
        """ANDOR has position_file_pattern=null and position_format=ims."""
        cfg = load_microscope_config("andor")
        assert cfg.position_file_pattern is None
        assert cfg.position_format == "ims"
        assert cfg.stage_x_heading is None
        assert cfg.stage_y_heading is None

    def test_invalid_microscope_name_raises(self):
        with pytest.raises(ValueError, match="Unknown microscope"):
            load_microscope_config("unknown_scope")
