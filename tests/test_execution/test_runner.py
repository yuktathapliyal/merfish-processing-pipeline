"""Tests for the execution runner — enabled flags and stage registry."""

import pytest

from merfish_pipeline.config.defaults import VALID_STAGES
from merfish_pipeline.stages.registry import STAGE_REGISTRY, list_stages


class TestStageRegistry:
    """Verify that all stages are registered."""

    def test_registry_has_all_stages(self):
        """All VALID_STAGES should have a registered implementation."""
        registered = set(STAGE_REGISTRY.keys())
        expected = set(VALID_STAGES)
        assert registered == expected, (
            f"Missing registrations: {expected - registered}; "
            f"Extra registrations: {registered - expected}"
        )

    def test_registry_count(self):
        """Registry count should match VALID_STAGES."""
        assert len(STAGE_REGISTRY) == len(VALID_STAGES)

    def test_list_stages_returns_sorted(self):
        """list_stages() should return all names sorted alphabetically."""
        names = list_stages()
        assert names == sorted(names)
        assert len(names) == len(VALID_STAGES)


class TestEnabledFlagSkipping:
    """Verify that stages with enabled=False are skipped by the runner."""

    def test_runner_checks_enabled_flag(self):
        """The runner loop should check for stage_cfg.enabled and skip if False.

        This is a structural test — we verify that the check exists in the
        runner source code.
        """
        from merfish_pipeline.execution import runner
        import inspect

        source = inspect.getsource(runner.run_pipeline)
        # The runner should check for an 'enabled' attribute
        assert "enabled" in source, (
            "run_pipeline should check stage_cfg.enabled"
        )
        # And it should skip disabled stages
        assert "skipped" in source.lower() or "skip" in source.lower(), (
            "run_pipeline should skip disabled stages"
        )
