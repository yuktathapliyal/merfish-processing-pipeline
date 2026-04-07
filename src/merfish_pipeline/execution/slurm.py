"""SLURM job generation for the merFISH processing pipeline.

Generates a master submission script that chains pipeline stages as
dependent ``sbatch`` jobs.  Each job invokes ``merfish-pipe run`` with
a single ``--stage`` flag, so stages execute serially on the cluster
with automatic dependency tracking.
"""

from __future__ import annotations

import logging
import shlex
from pathlib import Path
from typing import Any

logger = logging.getLogger("merfish_pipeline.execution.slurm")


def generate_slurm_script(
    config: Any,
    stages: list[str],
    experiment_yaml: Path,
    output_path: Path | None = None,
) -> Path:
    """Generate a SLURM submission script that chains pipeline stages.

    Each stage becomes its own ``sbatch`` job, chained via
    ``--dependency=afterok:$PREV_JOB_ID``.

    Parameters
    ----------
    config:
        Fully resolved ``PipelineConfig``.
    stages:
        Ordered list of stage names to submit.
    experiment_yaml:
        Path to the experiment config file (passed to ``merfish-pipe run``).
    output_path:
        Where to write the script.  Defaults to
        ``{output_dir}/submit_{experiment_name}.sh``.

    Returns
    -------
    Path
        The path to the generated script.
    """
    slurm_cfg = config.execution.slurm
    if slurm_cfg is None:
        raise ValueError(
            "SLURM config is not set.  Use --profile slurm or configure "
            "execution.slurm in your profile YAML."
        )

    xp_name = config.experiment.name
    logs_dir = config.logs_dir or Path(config.paths.output_dir) / "logs"
    experiment_yaml = Path(experiment_yaml).resolve()

    if output_path is None:
        output_path = Path(config.paths.output_dir) / f"submit_{xp_name}.sh"

    lines: list[str] = [
        "#!/bin/bash",
        f"# Auto-generated SLURM submission script for experiment: {xp_name}",
        f"# Stages: {', '.join(stages)}",
        "",
        f"LOG_DIR=\"{logs_dir}\"",
        "mkdir -p \"$LOG_DIR\"",
        "",
        "PREV_JOB_ID=\"\"",
        "",
    ]

    for i, stage_name in enumerate(stages):
        job_name = f"merfish_{xp_name}_{stage_name}"
        log_file = f"$LOG_DIR/slurm_{stage_name}_%j.log"

        sbatch_lines = [
            f"# --- Stage {i + 1}: {stage_name} ---",
            "SBATCH_ARGS=()",
            f'SBATCH_ARGS+=(--job-name="{job_name}")',
            f'SBATCH_ARGS+=(--output="{log_file}")',
            f'SBATCH_ARGS+=(--error="{log_file}")',
            f'SBATCH_ARGS+=(--partition="{slurm_cfg.partition}")',
            f'SBATCH_ARGS+=(--time="{slurm_cfg.time}")',
            f'SBATCH_ARGS+=(--mem-per-cpu="{slurm_cfg.mem_per_cpu}")',
            f'SBATCH_ARGS+=(--cpus-per-task={slurm_cfg.cpus_per_task})',
        ]

        # Dependency on previous job
        if i > 0:
            sbatch_lines.append(
                'if [ -n "$PREV_JOB_ID" ]; then'
            )
            sbatch_lines.append(
                '    SBATCH_ARGS+=(--dependency=afterok:$PREV_JOB_ID)'
            )
            sbatch_lines.append("fi")

        # The actual submission command.  ``shlex.quote`` ensures that paths
        # or stage names containing spaces / quotes / shell metacharacters
        # are passed safely to ``sbatch``.  Outer ``--wrap=`` quoting stays
        # double-quoted so the surrounding ``${SBATCH_ARGS[@]}`` expansion
        # and the ``awk`` subshell continue to work.
        yaml_quoted = shlex.quote(str(experiment_yaml))
        stage_quoted = shlex.quote(stage_name)
        cmd = (
            f"merfish-pipe run {yaml_quoted} --profile slurm "
            f"--stage {stage_quoted} --slurm-worker"
        )
        sbatch_lines.append(
            f'PREV_JOB_ID=$(sbatch "${{SBATCH_ARGS[@]}}" --wrap="{cmd}" | awk \'{{print $4}}\')'
        )
        sbatch_lines.append(f'echo "Submitted {stage_name}: job $PREV_JOB_ID"')
        sbatch_lines.append("")

        lines.extend(sbatch_lines)

    lines.append("echo \"\"")
    lines.append(f'echo "All {len(stages)} stages submitted for {xp_name}."')
    lines.append('echo "Monitor with: squeue -u $USER"')

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    output_path.chmod(0o755)

    logger.info("SLURM submission script written to: %s", output_path)
    return output_path
