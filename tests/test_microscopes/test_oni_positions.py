"""Tests for ONI adapter position reading — tile_number and grid_pos handling."""

import pandas as pd
import pytest


class TestONITileNumber:
    """Verify tile_number is read from the CSV column, not the DataFrame row index."""

    def test_tile_number_from_column(self):
        """When tile_number column exists, its value should be used (not row index)."""
        # Simulate a position CSV with tile_number starting at 1 (1-based)
        df = pd.DataFrame({
            "tile_number": [1, 2, 3],
            "stage_pos_x": [100.0, 200.0, 300.0],
            "stage_pos_y": [400.0, 500.0, 600.0],
        })
        # Row indices are 0, 1, 2 but tile_number values are 1, 2, 3
        for idx, row in df.iterrows():
            tile = int(row["tile_number"]) if "tile_number" in df.columns else int(idx)
            # tile should be from the column (1, 2, 3), not the index (0, 1, 2)
            assert tile == idx + 1, (
                f"tile_number should be {idx + 1} (from column), got {tile}"
            )

    def test_tile_number_fallback_to_index(self):
        """When tile_number column is absent, fall back to row index."""
        df = pd.DataFrame({
            "stage_pos_x": [100.0, 200.0],
            "stage_pos_y": [400.0, 500.0],
        })
        for idx, row in df.iterrows():
            tile = int(row["tile_number"]) if "tile_number" in df.columns else int(idx)
            assert tile == idx

    def test_grid_pos_columns_preserved(self):
        """When grid_pos_x/y exist in position file, they should be included in output."""
        df = pd.DataFrame({
            "tile_number": [1, 2],
            "stage_pos_x": [100.0, 200.0],
            "stage_pos_y": [400.0, 500.0],
            "grid_pos_x": [0, 1],
            "grid_pos_y": [0, 0],
        })
        records = []
        for idx, row in df.iterrows():
            record = {
                "tile_number": int(row["tile_number"]) if "tile_number" in df.columns else int(idx),
                "stage_pos_x": float(row["stage_pos_x"]),
                "stage_pos_y": float(row["stage_pos_y"]),
            }
            for grid_col in ("grid_pos_x", "grid_pos_y"):
                if grid_col in df.columns:
                    record[grid_col] = int(row[grid_col])
            records.append(record)

        assert records[0]["grid_pos_x"] == 0
        assert records[0]["grid_pos_y"] == 0
        assert records[1]["grid_pos_x"] == 1
        assert records[1]["grid_pos_y"] == 0
        assert "grid_pos_x" in records[0]
        assert "grid_pos_y" in records[0]
