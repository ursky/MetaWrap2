
"""Tests for global CLI flag propagation from the run/cohort driver to child modules.

Guards a past bug: the driver applied global resource and execution overrides to the
parent process but did not forward them to child module invocations. As a result,
modules such as read_qc used default limits (e.g. max_disk_fraction = 0.7) even when
the parent was launched with --max-disk-fraction 0.9.

The tests check resource overrides, execution flags, and the generated child commands.
"""

import pytest

from metawrap2 import command, limits
from metawrap2.commands import run as driver


@pytest.fixture(autouse=True)
def reset_global_flags():
    # Each test starts without overrides inherited from another test.
    limits.clear_overrides()
    command.set_skip_space_check(False)
    command.set_skip_validation(False)
    command.set_verbose_logs(False)
    command.set_heartbeat_seconds(None)

    yield

    limits.clear_overrides()
    command.set_skip_space_check(False)
    command.set_skip_validation(False)
    command.set_verbose_logs(False)
    command.set_heartbeat_seconds(None)


def _child_flags(force=False, resume=False, dry_run=False):
    return driver.child_global_flags(
        force=force,
        resume=resume,
        dry_run=dry_run,
    )

# --- resource limit overrides -------------------------------------------------------------


@pytest.mark.parametrize(
    "key,value,expected",
    [
        ("max_threads", 32, ["--max-threads", "32"]),
        ("max_thread_fraction", 0.8, ["--max-thread-fraction", "0.8"]),
        ("max_memory_gb", 1000, ["--max-memory-gb", "1000.0"]),
        ("max_memory_fraction", 0.9, ["--max-memory-fraction", "0.9"]),
        ("max_disk_fraction", 0.9, ["--max-disk-fraction", "0.9"]),
        ("max_gpus", 0, ["--max-gpus", "0"]),
        ("disk_check_interval", 30, ["--disk-check-interval", "30"]),
        ("gpu_strict", True, ["--gpu-strict"]),
    ],
)
def test_resource_overrides_are_forwarded(key, value, expected):
    limits.set_override(key, value)

    flags = _child_flags()

    assert flags == expected


def test_multiple_resource_overrides_are_forwarded_together():
    limits.set_override("max_threads", 32)
    limits.set_override("max_memory_gb", 1000)
    limits.set_override("max_disk_fraction", 0.9)

    flags = _child_flags()

    assert flags == [
        "--max-threads", "32",
        "--max-memory-gb", "1000.0",
        "--max-disk-fraction", "0.9",
    ]


def test_gpu_strict_false_is_not_forwarded():
    # --gpu-strict is a switch, not a flag that accepts a boolean value.
    limits.set_override("gpu_strict", False)

    assert "--gpu-strict" not in _child_flags()


# --- global execution flags ---------------------------------------------------------------


def test_execution_flags_are_forwarded():
    command.set_skip_space_check(True)
    command.set_skip_validation(True)
    command.set_verbose_logs(True)
    command.set_heartbeat_seconds(15)

    flags = _child_flags(
        force=True,
        resume=True,
        dry_run=True,
    )

    assert flags == [
        "--force",
        "--resume",
        "--dry-run",
        "--skip-space-check",
        "--skip-validation",
        "--verbose-logs",
        "--heartbeat-seconds", "15.0",
    ]


def test_unset_execution_flags_are_not_forwarded():
    # Child modules should not receive flags that the parent did not enable.
    assert _child_flags() == []


def test_heartbeat_zero_is_forwarded():
    # Zero explicitly disables heartbeats and must not be treated as unset.
    command.set_heartbeat_seconds(0)

    assert _child_flags() == [
        "--heartbeat-seconds", "0.0",
    ]


# --- child module commands ----------------------------------------------------------------


def test_child_command_includes_global_flags():
    limits.set_override("max_disk_fraction", 0.9)
    command.set_skip_space_check(True)

    flags = _child_flags()

    step = driver.Step(
        name="read_qc:A",
        module="read_qc",
        argv=[
            "-1", "raw/A_1.fastq.gz",
            "-2", "raw/A_2.fastq.gz",
            "-o", "study/READ_QC/A",
        ],
    )

    argv = step.command(
        threads=16,
        config=None,
        global_flags=flags,
    )

    assert "--max-disk-fraction" in argv
    assert argv[argv.index("--max-disk-fraction") + 1] == "0.9"
    assert "--skip-space-check" in argv
    assert "read_qc" in argv

    # Global flags must precede the module name for the CLI parser.
    assert argv.index("--max-disk-fraction") < argv.index("read_qc")
    assert argv.index("--skip-space-check") < argv.index("read_qc")

    # Per-step thread allocation must still be passed to the child.
    assert argv[argv.index("-t") + 1] == "16"


def test_child_command_preserves_module_arguments():
    flags = _child_flags(resume=True)

    step = driver.Step(
        name="assembly:A",
        module="assembly",
        argv=[
            "-1", "reads/A_1.fastq.gz",
            "-2", "reads/A_2.fastq.gz",
            "-o", "study/ASSEMBLY/A",
        ],
    )

    argv = step.command(
        threads=16,
        config="metawrap2.toml",
        global_flags=flags,
        memory_gb=60,
    )

    assert "--resume" in argv
    assert "assembly" in argv
    assert argv[argv.index("-t") + 1] == "16"
    assert argv[argv.index("-m", argv.index("assembly")) + 1] == "60"
    assert argv[argv.index("--config") + 1] == "metawrap2.toml"



# --- workflow entry points ----------------------------------------------------------------


@pytest.mark.parametrize("cohort", [False, True])
def test_workflow_forwards_global_flags(tmp_path, monkeypatch, cohort):
    """Both run and cohort must pass global overrides to the workflow executor."""
    from metawrap2.commands import cohort as cohort_cmd

    sheet = tmp_path / "samples.toml"
    sheet.write_text("""
[[samples]]
name = "A"
r1 = "raw/A_1.fastq.gz"
r2 = "raw/A_2.fastq.gz"
""")

    seen = {}

    def fake_execute(
        steps, layout, study, threads, jobs, config,
        global_flags, state, dry_run, **kwargs
    ):
        seen["flags"] = global_flags
        seen["steps"] = steps
        return 0

    monkeypatch.setattr(driver, "execute", fake_execute)

    # Simulate global CLI overrides applied by metawrap2.cli.
    limits.set_override("max_disk_fraction", 0.9)
    command.set_skip_space_check(True)
    command.set_skip_validation(True)

    argv = [
        str(sheet),
        "-o", str(tmp_path / "study"),
    ]

    if cohort:
        result = cohort_cmd.main(argv)
    else:
        result = driver.main(argv)

    assert result == 0

    flags = seen["flags"]

    assert "--max-disk-fraction" in flags
    assert flags[flags.index("--max-disk-fraction") + 1] == "0.9"
    assert "--skip-space-check" in flags
    assert "--skip-validation" in flags

    # Verify that the generated flags reach the child module command.
    read_qc = next(
        step for step in seen["steps"]
        if step.module == "read_qc"
    )

    child_argv = read_qc.command(
        threads=16,
        config=None,
        global_flags=flags,
    )

    assert "--max-disk-fraction" in child_argv
    assert "--skip-space-check" in child_argv
    assert "--skip-validation" in child_argv

    # Global flags must precede the module name.
    module_index = child_argv.index("read_qc")

    assert child_argv.index("--max-disk-fraction") < module_index
    assert child_argv.index("--skip-space-check") < module_index
    assert child_argv.index("--skip-validation") < module_index
