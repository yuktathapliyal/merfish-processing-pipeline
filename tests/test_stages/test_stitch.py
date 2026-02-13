"""Tests for stitch stage — _build_position_map grid computation."""

import pandas as pd
import pytest

from merfish_pipeline.stages.stitch import _build_position_map


class TestBuildPositionMap:
    """Verify _build_position_map produces correct grid indices."""

    def test_grid_pos_columns_used_directly(self):
        """When grid_pos_x/y exist, they should be used as grid coords."""
        df = pd.DataFrame({
            "tile_number": [1, 2, 3, 4],
            "grid_pos_x": [0, 1, 0, 1],
            "grid_pos_y": [0, 0, 1, 1],
            "stage_pos_x": [0.0, 100.0, 0.0, 100.0],
            "stage_pos_y": [0.0, 0.0, 100.0, 100.0],
        })
        pos_map = _build_position_map(df)

        assert pos_map[1] == (0, 0)
        assert pos_map[2] == (1, 0)
        assert pos_map[3] == (0, 1)
        assert pos_map[4] == (1, 1)

    def test_grid_pos_shifted_to_zero(self):
        """Grid positions should be shifted so minimum becomes 0."""
        df = pd.DataFrame({
            "tile_number": [1, 2],
            "grid_pos_x": [5, 6],
            "grid_pos_y": [10, 10],
        })
        pos_map = _build_position_map(df)

        assert pos_map[1] == (0, 0)
        assert pos_map[2] == (1, 0)

    def test_stage_pos_fallback_produces_small_grid(self):
        """When grid_pos columns are absent, stage positions should be
        quantized to small grid indices (not raw micron values)."""
        # Simulating a 3x2 grid with 100um pitch
        df = pd.DataFrame({
            "tile_number": [1, 2, 3, 4, 5, 6],
            "stage_pos_x": [0.0, 100.0, 200.0, 0.0, 100.0, 200.0],
            "stage_pos_y": [0.0, 0.0, 0.0, 100.0, 100.0, 100.0],
        })
        pos_map = _build_position_map(df)

        # Grid indices should be small integers (0, 1, 2), not 0, 100, 200
        max_gx = max(gx for gx, _ in pos_map.values())
        max_gy = max(gy for _, gy in pos_map.values())

        assert max_gx == 2, f"Expected max grid_x=2, got {max_gx}"
        assert max_gy == 1, f"Expected max grid_y=1, got {max_gy}"

        # Specific positions
        assert pos_map[1] == (0, 0)
        assert pos_map[2] == (1, 0)
        assert pos_map[3] == (2, 0)
        assert pos_map[4] == (0, 1)
        assert pos_map[5] == (1, 1)
        assert pos_map[6] == (2, 1)

    def test_stage_pos_fallback_single_tile(self):
        """Single tile should get grid position (0, 0)."""
        df = pd.DataFrame({
            "tile_number": [1],
            "stage_pos_x": [500.0],
            "stage_pos_y": [300.0],
        })
        pos_map = _build_position_map(df)
        assert pos_map[1] == (0, 0)

    def test_stage_pos_fallback_irregular_spacing(self):
        """With irregular spacing, the minimum gap should be used as pitch."""
        # tiles at x=0, 50, 100, 200 — min gap is 50
        df = pd.DataFrame({
            "tile_number": [1, 2, 3, 4],
            "stage_pos_x": [0.0, 50.0, 100.0, 200.0],
            "stage_pos_y": [0.0, 0.0, 0.0, 0.0],
        })
        pos_map = _build_position_map(df)

        assert pos_map[1] == (0, 0)
        assert pos_map[2] == (1, 0)
        assert pos_map[3] == (2, 0)
        assert pos_map[4] == (4, 0)
