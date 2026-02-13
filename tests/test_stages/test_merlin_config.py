"""Tests for merlin_config stage — positions format and run script generation."""

import json
from pathlib import Path

import pandas as pd
import pytest

from merfish_pipeline.stages.merlin_config import (
    _generate_positions_csv,
    _generate_run_script,
)


class TestGeneratePositionsCsv:
    """Verify MERlin positions CSV is 2-column, headerless."""

    def test_two_column_no_header(self, tmp_path):
        """Output should have exactly 2 columns (x, y) and no header row."""
        positions_df = pd.DataFrame({
            "round": [1, 1, 1],
            "tile_number": [1, 2, 3],
            "stage_pos_x": [100.5, 200.3, 300.1],
            "stage_pos_y": [400.2, 500.4, 600.6],
        })
        out_path = tmp_path / "positions.csv"
        _generate_positions_csv(positions_df, out_path)

        text = out_path.read_text()
        lines = [l for l in text.strip().split("\n") if l]

        assert len(lines) == 3, f"Expected 3 data lines, got {len(lines)}"

        # First line should NOT be a header — it should be numeric
        first_parts = lines[0].split(",")
        assert len(first_parts) == 2, f"Expected 2 columns, got {len(first_parts)}"
        # Should be parseable as floats (not column names)
        float(first_parts[0])
        float(first_parts[1])

    def test_sorted_by_tile_number(self, tmp_path):
        """Positions should be ordered by tile_number."""
        positions_df = pd.DataFrame({
            "round": [1, 1, 1],
            "tile_number": [3, 1, 2],
            "stage_pos_x": [300.0, 100.0, 200.0],
            "stage_pos_y": [600.0, 400.0, 500.0],
        })
        out_path = tmp_path / "positions.csv"
        _generate_positions_csv(positions_df, out_path)

        text = out_path.read_text()
        lines = [l for l in text.strip().split("\n") if l]
        x_values = [float(l.split(",")[0]) for l in lines]
        assert x_values == [100.0, 200.0, 300.0]

    def test_empty_df_writes_empty_file(self, tmp_path):
        """Empty input should produce an empty file."""
        out_path = tmp_path / "positions.csv"
        _generate_positions_csv(pd.DataFrame(), out_path)
        assert out_path.read_text() == ""


class TestGenerateRunScript:
    """Verify the run_merLIN.sh script has correct MERlin command syntax."""

    def test_uses_o_flag_not_d(self, tmp_path):
        """Data organisation should use -o (not -d)."""
        merlin_dir = tmp_path / "merlin_data"
        merlin_dir.mkdir()
        out = merlin_dir / "run_merLIN.sh"

        _generate_run_script(
            merlin_data_dir=merlin_dir,
            xp_name="XP001",
            codebook_filename="codebook.csv",
            cores=100,
            output_path=out,
        )
        script = out.read_text()

        assert "-o data_organization_XP001.csv" in script, (
            "Script should use -o for data organisation"
        )
        assert "-d " not in script, "Script should NOT use -d flag"

    def test_has_merlin_env_path_export(self, tmp_path):
        """Script should export MERLIN_ENV_PATH."""
        merlin_dir = tmp_path / "merlin_data"
        merlin_dir.mkdir()
        out = merlin_dir / "run_merLIN.sh"

        _generate_run_script(
            merlin_data_dir=merlin_dir,
            xp_name="XP001",
            codebook_filename="codebook.csv",
            cores=100,
            output_path=out,
        )
        script = out.read_text()

        assert "export MERLIN_ENV_PATH=" in script, (
            "Script should set MERLIN_ENV_PATH"
        )

    def test_has_trailing_data_dir(self, tmp_path):
        """MERlin command should end with the data directory argument."""
        merlin_dir = tmp_path / "merlin_data"
        merlin_dir.mkdir()
        out = merlin_dir / "run_merLIN.sh"

        _generate_run_script(
            merlin_data_dir=merlin_dir,
            xp_name="XP001",
            codebook_filename="codebook.csv",
            cores=100,
            output_path=out,
        )
        script = out.read_text()

        # The trailing argument should be the directory name
        assert "merlin_data" in script, (
            "Script should include trailing data directory argument"
        )

    def test_no_source_activate(self, tmp_path):
        """Script should NOT contain 'source activate' — env activation is user's job."""
        merlin_dir = tmp_path / "merlin_data"
        merlin_dir.mkdir()
        out = merlin_dir / "run_merLIN.sh"

        _generate_run_script(
            merlin_data_dir=merlin_dir,
            xp_name="XP001",
            codebook_filename="codebook.csv",
            cores=100,
            output_path=out,
        )
        script = out.read_text()

        assert "source activate" not in script

    def test_script_is_executable(self, tmp_path):
        """Output script should have executable permissions."""
        merlin_dir = tmp_path / "merlin_data"
        merlin_dir.mkdir()
        out = merlin_dir / "run_merLIN.sh"

        _generate_run_script(
            merlin_data_dir=merlin_dir,
            xp_name="XP001",
            codebook_filename="codebook.csv",
            cores=100,
            output_path=out,
        )

        import os
        assert os.access(out, os.X_OK)


class TestReregistrationDetection:
    """Verify merlin_config auto-detects reregistration output."""

    def test_detect_target_z_from_metadata(self, tmp_path):
        """When reregistration metadata exists, _detect_reregistration should
        return target_z and remapped_data_dir."""
        from merfish_pipeline.stages.merlin_config import MerlinConfigStage

        # Create a fake reregistration run_metadata.json
        rereg_dir = tmp_path / "output" / "reregistration"
        rereg_dir.mkdir(parents=True)
        metadata = {
            "parameters": {
                "total_z": 15,
                "target_z": 11,
                "min_start": 5,
            }
        }
        (rereg_dir / "run_metadata.json").write_text(
            json.dumps(metadata), encoding="utf-8"
        )

        # Build a minimal mock config
        class _MockRereg:
            enabled = True

        class _MockPaths:
            output_dir = tmp_path / "output"
            remapped_data_dir = tmp_path / "remapped_data"

        class _MockConfig:
            reregistration = _MockRereg()
            paths = _MockPaths()

        # Instantiate stage with mock (bypass __init__ via __new__)
        stage = object.__new__(MerlinConfigStage)
        stage.config = _MockConfig()

        import logging
        stage.logger = logging.getLogger("test")

        result = stage._detect_reregistration()
        assert result is not None
        assert result["target_z"] == 11
        assert result["remapped_data_dir"] == tmp_path / "remapped_data"

    def test_detect_returns_none_when_disabled(self, tmp_path):
        """When reregistration is disabled, should return None."""
        from merfish_pipeline.stages.merlin_config import MerlinConfigStage

        class _MockRereg:
            enabled = False

        class _MockConfig:
            reregistration = _MockRereg()

        stage = object.__new__(MerlinConfigStage)
        stage.config = _MockConfig()

        import logging
        stage.logger = logging.getLogger("test")

        assert stage._detect_reregistration() is None

    def test_detect_returns_none_when_no_metadata(self, tmp_path):
        """When enabled but no metadata file exists, should return None."""
        from merfish_pipeline.stages.merlin_config import MerlinConfigStage

        class _MockRereg:
            enabled = True

        class _MockPaths:
            output_dir = tmp_path / "output"  # no reregistration subdir

        class _MockConfig:
            reregistration = _MockRereg()
            paths = _MockPaths()

        stage = object.__new__(MerlinConfigStage)
        stage.config = _MockConfig()

        import logging
        stage.logger = logging.getLogger("test")

        assert stage._detect_reregistration() is None
