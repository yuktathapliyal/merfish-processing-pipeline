"""``merlin_config`` stage -- generate all MERlin parameter files.

This stage generates the configuration files required to run MERlin, but does
**not** invoke MERlin itself.  The outputs include:

1. **data_organization CSV** -- Frame-to-bit mapping expanded from a template.
2. **microscope JSON** -- Microscope parameters copied from a template.
3. **positions CSV** -- FOV tile positions in MERlin format.
4. **analysis JSON** -- Analysis parameters copied from a template.
5. **codebook CSV** -- Codebook copied from a template.
6. **.merlinenv** -- Environment file pointing to data/analysis/parameters dirs.
7. **run_merLIN.sh** -- Shell script to invoke MERlin.
"""

from __future__ import annotations

import json
import shutil
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from merfish_pipeline.io.sheet_io import read_sheet, write_sheet
from merfish_pipeline.stages.base import PipelineStage, StageResult
from merfish_pipeline.stages.registry import register_stage

# ---------------------------------------------------------------------------
# Package-level templates directory (two levels up from the ``src/`` tree)
# ---------------------------------------------------------------------------

_PACKAGE_DIR = Path(__file__).resolve().parent  # .../stages/
_TEMPLATES_DIR = _PACKAGE_DIR.parents[2] / "templates"  # repo-root/templates


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _resolve_template(path: Path | str, templates_dir: Path | None = None) -> Path:
    """Locate a template file.

    Resolution order:

    1. If *path* is absolute and exists, return it directly.
    2. If *templates_dir* is given, look for ``templates_dir / path`` and common
       sub-directories (``dataorganization/``, ``microscope/``, ``analysis/``,
       ``codebooks/``).
    3. Fall back to the package-level ``templates/`` directory using the same
       sub-directory search.

    Raises
    ------
    FileNotFoundError
        If the template cannot be found in any of the searched locations.
    """
    path = Path(path)

    # 1. Absolute and exists
    if path.is_absolute() and path.exists():
        return path

    search_dirs: list[Path] = []
    for base in [templates_dir, _TEMPLATES_DIR]:
        if base is None or not base.is_dir():
            continue
        search_dirs.append(base)
        for subdir in ("dataorganization", "microscope", "analysis", "codebooks"):
            candidate = base / subdir
            if candidate.is_dir():
                search_dirs.append(candidate)

    for directory in search_dirs:
        candidate = directory / path.name
        if candidate.exists():
            return candidate
        # Also try the path as-is (may include subdirectory components)
        candidate = directory / path
        if candidate.exists():
            return candidate

    searched = ", ".join(str(d) for d in search_dirs) if search_dirs else "(none)"
    raise FileNotFoundError(
        f"Template not found: {path}. Searched in: {searched}"
    )


def _extract_base_frame(value: object) -> int:
    """Extract 1-based frame from an expanded array string like ``[2 5 8 ...]``.

    The first element of the array equals ``original_frame - 1``, so we add 1
    to recover the original 1-based value.
    """
    inner = str(value).strip().strip("[]").split()
    return int(inner[0]) + 1


def _extract_base_fiducial(value: object) -> int:
    """Convert a 0-based fiducialFrame back to 1-based."""
    return int(value) + 1


def _expand_data_organization(
    template_df: pd.DataFrame,
    n_z: int,
    n_channels: int,
) -> pd.DataFrame:
    """Expand a data-organization template into the full MERlin format.

    The template CSV has one row per bit/readout.  The ``frame`` and
    ``fiducialFrame`` columns may contain either:

    * **Scalar** 1-based integers (non-expanded template), or
    * **Array strings** like ``[2 5 8 ...]`` (pre-expanded for a different
      experiment's z-count).

    In both cases this function (re-)expands them for the given *n_z* and
    *n_channels*:

    * ``frame`` becomes a numpy-style array string:
      ``[f0 f1 f2 ...]`` where ``f_z = (original_frame - 1) + z * n_channels``
      for each *z* in ``range(n_z)``.
    * ``zPos`` becomes ``[0. 1. 2. ... (n_z - 1).]``
    * ``fiducialFrame`` is converted from 1-indexed to 0-indexed (single int).

    Parameters
    ----------
    template_df:
        DataFrame read from the data-organization template CSV.
    n_z:
        Number of z-slices per FOV.
    n_channels:
        Number of colour channels per z-slice.

    Returns
    -------
    pd.DataFrame
        Expanded data organisation ready for MERlin.
    """
    df = template_df.copy()

    # If the template is pre-expanded (array-style frame values), normalise
    # back to scalar 1-based integers so the expansion below works correctly.
    first_frame = str(df["frame"].iloc[0]).strip()
    if "[" in first_frame:
        df["frame"] = df["frame"].apply(_extract_base_frame)
        df["fiducialFrame"] = df["fiducialFrame"].apply(_extract_base_fiducial)

    # Build the zPos array string (shared by all rows).
    z_positions = np.arange(n_z, dtype=float)
    z_pos_str = np.array2string(z_positions, separator=" ", max_line_width=10000)

    expanded_frames: list[str] = []
    expanded_fid_frames: list[int] = []

    for _, row in df.iterrows():
        original_frame = int(row["frame"])
        # frame = (original_frame - 1) + z * n_channels for each z
        frames = [(original_frame - 1) + z * n_channels for z in range(n_z)]
        frame_arr = np.array(frames, dtype=int)
        frame_str = np.array2string(frame_arr, separator=" ", max_line_width=10000)
        expanded_frames.append(frame_str)

        original_fid = int(row["fiducialFrame"])
        expanded_fid_frames.append(original_fid - 1)

    df["frame"] = expanded_frames
    df["zPos"] = z_pos_str
    df["fiducialFrame"] = expanded_fid_frames

    return df


def _generate_positions_csv(positions_df: pd.DataFrame, output_path: Path) -> None:
    """Write FOV tile positions in the MERlin-expected format.

    MERlin expects a headerless, two-column CSV (x, y) with one row per FOV,
    ordered by tile number.  The input *positions_df* comes from the ``index``
    stage with columns ``round, tile_number, stage_pos_x, stage_pos_y, ...``.

    We take one representative position per tile (from the first round
    encountered) and write only ``stage_pos_x, stage_pos_y``.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    if positions_df.empty:
        output_path.write_text("", encoding="utf-8")
        return

    # Take positions from the first round to avoid duplicates.
    if "round" in positions_df.columns:
        first_round = positions_df["round"].min()
        pos = positions_df[positions_df["round"] == first_round].copy()
    else:
        pos = positions_df.copy()

    # Determine the tile/fov column for sorting.
    if "tile_number" in pos.columns:
        fov_col = "tile_number"
    elif "fov" in pos.columns:
        fov_col = "fov"
    else:
        pos = pos.reset_index(drop=True)
        pos["fov"] = pos.index
        fov_col = "fov"

    pos = pos.sort_values(fov_col).reset_index(drop=True)

    result = pd.DataFrame(
        {
            "x": pos["stage_pos_x"].values,
            "y": pos["stage_pos_y"].values,
        }
    )
    result.to_csv(output_path, index=False, header=False)


def _generate_merlinenv(
    data_home: Path,
    analysis_home: Path,
    parameters_home: Path,
    output_path: Path,
) -> None:
    """Write a ``.merlinenv`` file.

    Format::

        DATA_HOME=/path/to/merlin_data
        ANALYSIS_HOME=/path/to/merlin_analysis
        PARAMETERS_HOME=/path/to/parameters
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"DATA_HOME={data_home}",
        f"ANALYSIS_HOME={analysis_home}",
        f"PARAMETERS_HOME={parameters_home}",
    ]
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _generate_run_script(
    merlin_data_dir: Path,
    merlinenv_dir: Path,
    xp_name: str,
    codebook_filename: str,
    cores: int,
    output_path: Path,
) -> None:
    """Write a ``run_merLIN.sh`` shell script that invokes MERlin.

    Parameters
    ----------
    merlin_data_dir:
        Directory containing the merged TIFF data.  Its name is passed as
        the positional argument to MERlin.
    merlinenv_dir:
        Directory containing the ``.merlinenv`` file.  ``MERLIN_ENV_PATH``
        is set to this directory.
    xp_name:
        Experiment name used in parameter file naming.
    codebook_filename:
        Filename of the codebook CSV (not the full path).
    cores:
        Number of cores to pass to MERlin's ``-n`` flag.
    output_path:
        Where to write the shell script.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # The trailing positional argument is the data directory name relative to
    # DATA_HOME.
    data_dir_rel = merlin_data_dir.name

    script = f"""\
#!/bin/bash
export MERLIN_ENV_PATH="{merlinenv_dir}"

merlin \\
    -a analysis_{xp_name}.json \\
    -c {codebook_filename} \\
    -o data_organization_{xp_name}.csv \\
    -m microscope_{xp_name}.json \\
    -p positions_{xp_name}.csv \\
    -n {cores} \\
    {data_dir_rel}
"""
    output_path.write_text(script, encoding="utf-8")
    # Make executable
    output_path.chmod(output_path.stat().st_mode | 0o755)


# ---------------------------------------------------------------------------
# Stage implementation
# ---------------------------------------------------------------------------


@register_stage("merlin_config")
class MerlinConfigStage(PipelineStage):
    """Generate all MERlin parameter files (does not run MERlin)."""

    description = "Generate MERlin configuration and parameter files"

    # ------------------------------------------------------------------
    # Abstract interface
    # ------------------------------------------------------------------

    def validate_inputs(self) -> list[str]:
        errors: list[str] = []

        # Manifest from index stage
        manifest_path = self._manifest_path()
        if not manifest_path.exists():
            errors.append(
                f"Manifest not found (run 'index' stage first): {manifest_path}"
            )

        # Positions from index stage
        positions_path = self._positions_path()
        if not positions_path.exists():
            errors.append(
                f"Positions file not found (run 'index' stage first): {positions_path}"
            )

        # Data organisation template
        try:
            self._resolve_data_org_template()
        except FileNotFoundError as exc:
            errors.append(str(exc))

        # Codebook template
        if self.config.merlin.codebook_template is not None:
            try:
                _resolve_template(self.config.merlin.codebook_template)
            except FileNotFoundError as exc:
                errors.append(str(exc))

        # Analysis template
        if self.config.merlin.analysis_template is not None:
            try:
                _resolve_template(self.config.merlin.analysis_template)
            except FileNotFoundError as exc:
                errors.append(str(exc))

        # Microscope template
        if self.config.merlin.microscope_template is not None:
            try:
                _resolve_template(self.config.merlin.microscope_template)
            except FileNotFoundError as exc:
                errors.append(str(exc))

        return errors

    def check_outputs_exist(self) -> bool:
        """Return True if all key output files already exist."""
        params_dir = Path(self.config.paths.parameters_dir)
        analysis_home = Path(self.config.paths.output_dir) / "merlin_analysis"
        xp_name = self.config.experiment.name

        required = [
            params_dir / "dataorganization" / f"data_organization_{xp_name}.csv",
            params_dir / "positions" / f"positions_{xp_name}.csv",
            analysis_home / ".merlinenv",
            analysis_home / "run_merLIN.sh",
            self.get_output_dir() / "run_metadata.json",
        ]
        return all(p.exists() for p in required)

    def run(self, dry_run: bool = False) -> StageResult:
        start_time = datetime.now()
        output_dir = self.get_output_dir()
        output_dir.mkdir(parents=True, exist_ok=True)

        xp_name = self.config.experiment.name
        params_dir = Path(self.config.paths.parameters_dir)
        merlin_data_dir = Path(self.config.paths.merlin_data_dir)
        analysis_home = Path(self.config.paths.output_dir) / "merlin_analysis"

        output_files: list[str] = []

        # ---------------------------------------------------------------
        # 0. Read manifest to get n_z_slices and n_channels
        # ---------------------------------------------------------------
        manifest_path = self._manifest_path()
        self.logger.info("Reading manifest from %s", manifest_path)
        manifest_df = read_sheet(manifest_path)

        n_z_raw = manifest_df["z_slice"].nunique()
        n_z = n_z_raw
        n_channels = manifest_df["channel"].nunique()
        self.logger.info(
            "Manifest: %d z-slices, %d channels", n_z_raw, n_channels,
        )

        # Check if reregistration was run — if so, override n_z and data_home.
        rereg_info = self._detect_reregistration()
        if rereg_info is not None:
            n_z = rereg_info["target_z"]
            self.logger.info(
                "Reregistration detected: using target_z=%d instead of raw n_z=%d",
                n_z,
                n_z_raw,
            )

        if dry_run:
            return StageResult(
                status="skipped",
                metadata={
                    "dry_run": True,
                    "n_z_slices": n_z,
                    "n_channels": n_channels,
                    "experiment_name": xp_name,
                },
            )

        # ---------------------------------------------------------------
        # 1. Data organisation CSV
        # ---------------------------------------------------------------
        data_org_template_path = self._resolve_data_org_template()
        self.logger.info(
            "Reading data organisation template: %s", data_org_template_path,
        )
        template_df = read_sheet(data_org_template_path)

        expanded_df = _expand_data_organization(template_df, n_z, n_channels)

        data_org_dir = params_dir / "dataorganization"
        data_org_dir.mkdir(parents=True, exist_ok=True)
        data_org_path = data_org_dir / f"data_organization_{xp_name}.csv"
        write_sheet(expanded_df, data_org_path)
        self.logger.info("Wrote data organisation: %s (%d rows)", data_org_path, len(expanded_df))
        output_files.append(str(data_org_path))

        # ---------------------------------------------------------------
        # 2. Microscope JSON
        # ---------------------------------------------------------------
        if self.config.merlin.microscope_template is not None:
            micro_src = _resolve_template(self.config.merlin.microscope_template)
            micro_dir = params_dir / "microscope"
            micro_dir.mkdir(parents=True, exist_ok=True)
            micro_dest = micro_dir / f"microscope_{xp_name}.json"
            shutil.copy2(micro_src, micro_dest)
            self.logger.info("Copied microscope parameters: %s -> %s", micro_src, micro_dest)
            output_files.append(str(micro_dest))

        # ---------------------------------------------------------------
        # 3. Positions CSV
        # ---------------------------------------------------------------
        positions_input_path = self._positions_path()
        self.logger.info("Reading positions from %s", positions_input_path)
        positions_df = read_sheet(positions_input_path)

        positions_dir = params_dir / "positions"
        positions_dir.mkdir(parents=True, exist_ok=True)
        positions_out = positions_dir / f"positions_{xp_name}.csv"
        _generate_positions_csv(positions_df, positions_out)
        self.logger.info("Wrote positions: %s", positions_out)
        output_files.append(str(positions_out))

        # ---------------------------------------------------------------
        # 4. Analysis JSON
        # ---------------------------------------------------------------
        analysis_filename: str | None = None
        if self.config.merlin.analysis_template is not None:
            analysis_src = _resolve_template(self.config.merlin.analysis_template)
            analysis_dir = params_dir / "analysis"
            analysis_dir.mkdir(parents=True, exist_ok=True)
            analysis_filename = f"analysis_{xp_name}.json"
            analysis_dest = analysis_dir / analysis_filename
            shutil.copy2(analysis_src, analysis_dest)
            self.logger.info("Copied analysis parameters: %s -> %s", analysis_src, analysis_dest)
            output_files.append(str(analysis_dest))

        # ---------------------------------------------------------------
        # 5. Codebook CSV
        # ---------------------------------------------------------------
        codebook_filename: str | None = None
        if self.config.merlin.codebook_template is not None:
            codebook_src = _resolve_template(self.config.merlin.codebook_template)
            codebook_dir = params_dir / "codebooks"
            codebook_dir.mkdir(parents=True, exist_ok=True)
            codebook_filename = codebook_src.name
            codebook_dest = codebook_dir / codebook_filename
            shutil.copy2(codebook_src, codebook_dest)
            self.logger.info("Copied codebook: %s -> %s", codebook_src, codebook_dest)
            output_files.append(str(codebook_dest))

        # ---------------------------------------------------------------
        # 6. .merlinenv
        # ---------------------------------------------------------------
        analysis_home.mkdir(parents=True, exist_ok=True)
        merlinenv_path = analysis_home / ".merlinenv"

        # DATA_HOME is the parent directory; MERlin resolves data at
        # DATA_HOME/<positional_arg> (e.g. output_dir/merlin_data).
        data_home = Path(self.config.paths.output_dir)

        if rereg_info is not None:
            data_dir = rereg_info["remapped_data_dir"]
            self.logger.info(
                "Using remapped data dir for MERlin: %s", data_dir,
            )
        else:
            data_dir = merlin_data_dir

        _generate_merlinenv(
            data_home=data_home,
            analysis_home=analysis_home,
            parameters_home=params_dir,
            output_path=merlinenv_path,
        )
        self.logger.info("Wrote .merlinenv: %s", merlinenv_path)
        output_files.append(str(merlinenv_path))

        # ---------------------------------------------------------------
        # 7. run_merLIN.sh
        # ---------------------------------------------------------------
        run_script_path = analysis_home / "run_merLIN.sh"
        _generate_run_script(
            merlin_data_dir=data_dir,
            merlinenv_dir=analysis_home,
            xp_name=xp_name,
            codebook_filename=codebook_filename or "codebook.csv",
            cores=self.config.merlin.cores,
            output_path=run_script_path,
        )
        self.logger.info("Wrote run script: %s", run_script_path)
        output_files.append(str(run_script_path))

        # ---------------------------------------------------------------
        # Summary and metadata
        # ---------------------------------------------------------------
        result = StageResult(
            status="completed",
            output_files=output_files,
            metadata={
                "experiment_name": xp_name,
                "n_z_slices": n_z,
                "n_z_raw": n_z_raw,
                "reregistration_applied": rereg_info is not None,
                "n_channels": n_channels,
                "data_org_rows": len(expanded_df),
                "n_positions": len(positions_df),
                "parameters_dir": str(params_dir),
                "merlin_data_dir": str(merlin_data_dir),
                "data_home": str(data_home),
            },
        )

        self.write_run_metadata(result, start_time)
        return result

    # ------------------------------------------------------------------
    # Internal: reregistration detection
    # ------------------------------------------------------------------

    def _detect_reregistration(self) -> dict | None:
        """Check if the reregistration stage ran and return its metadata.

        Returns
        -------
        dict or None
            If reregistration ran successfully, returns a dict with keys
            ``target_z`` (int) and ``remapped_data_dir`` (Path).
            Returns ``None`` if reregistration was not enabled or has no output.
        """
        if not self.config.reregistration.enabled:
            return None

        rereg_metadata_path = (
            Path(self.config.paths.output_dir) / "reregistration" / "run_metadata.json"
        )
        if not rereg_metadata_path.exists():
            return None

        try:
            with open(rereg_metadata_path, "r", encoding="utf-8") as fh:
                metadata = json.load(fh)
            target_z = metadata.get("parameters", {}).get("target_z")
            if target_z is None:
                self.logger.warning(
                    "Reregistration metadata exists but missing target_z: %s",
                    rereg_metadata_path,
                )
                return None
            remapped_data_dir = Path(self.config.paths.remapped_data_dir)
            return {"target_z": int(target_z), "remapped_data_dir": remapped_data_dir}
        except (json.JSONDecodeError, KeyError) as exc:
            self.logger.warning(
                "Could not parse reregistration metadata (%s): %s",
                rereg_metadata_path,
                exc,
            )
            return None

    # ------------------------------------------------------------------
    # Internal: path resolution helpers
    # ------------------------------------------------------------------

    def _manifest_path(self) -> Path:
        """Path to the manifest CSV produced by the index stage."""
        return Path(self.config.paths.output_dir) / "index" / "manifest.csv"

    def _positions_path(self) -> Path:
        """Path to the standardized positions CSV produced by the index stage."""
        return Path(self.config.paths.output_dir) / "index" / "positions.standardized.csv"

    def _resolve_data_org_template(self) -> Path:
        """Resolve the data-organisation template path.

        Uses ``config.raw_data.data_org_template`` and falls back to the
        package templates directory.
        """
        template_path = self.config.raw_data.data_org_template
        if template_path is None:
            raise FileNotFoundError(
                "No data_org_template configured in raw_data section."
            )
        return _resolve_template(template_path)
