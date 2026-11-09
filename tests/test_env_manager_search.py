"""Finding mamba or conda when neither is on PATH.

A non-login shell or a batch job often starts without conda on PATH, and the command runner then
could not start any tool in a conda environment. It now also looks in ``[settings] env_manager``,
``$MAMBA_EXE``, ``$CONDA_EXE`` and the usual install locations.
"""

from __future__ import annotations

import os

import pytest

from metawrap2 import command
from metawrap2.config import SETTINGS_KEYS, load_settings


def _fake_program(path):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("#!/bin/sh\nexit 0\n")
    path.chmod(0o755)
    return str(path)


@pytest.fixture
def bare_path(tmp_path, monkeypatch):
    """PATH with nothing on it, and the runner's manager restored afterwards."""
    empty = tmp_path / "empty_bin"
    empty.mkdir()
    monkeypatch.setenv("PATH", str(empty))
    monkeypatch.setattr(command.runner, "env_manager", "mamba")
    monkeypatch.setattr(command.runner, "dry_run", False)
    return tmp_path


def test_without_any_manager_the_message_names_what_to_set(bare_path):
    assert command.locate_env_manager("mamba") is None
    with pytest.raises(command.EnvManagerMissing) as caught:
        command.runner.env_prefix("metawrap2-binqc")
    text = str(caught.value)
    assert "MAMBA_EXE" in text and "[settings] env_manager" in text and "use_conda_envs" in text


def test_an_install_in_a_usual_location_is_found(bare_path, monkeypatch):
    root = bare_path / "miniforge3"
    mamba = _fake_program(root / "bin" / "mamba")
    monkeypatch.setattr(command, "ENV_MANAGER_ROOTS", (str(bare_path / "nothing"), str(root)))
    monkeypatch.setattr(command, "_RUN_FLAGS_CACHE", {mamba: []})
    assert command.runner.env_prefix("metawrap2-binqc") == [mamba, "run", "-n", "metawrap2-binqc"]
    assert command.runner.env_manager == mamba


def test_conda_is_used_when_there_is_no_mamba(bare_path, monkeypatch):
    conda = _fake_program(bare_path / "miniconda3" / "condabin" / "conda")
    monkeypatch.setattr(command, "ENV_MANAGER_ROOTS", (str(bare_path / "miniconda3"),))
    assert command.locate_env_manager("mamba") == conda


def test_the_environment_variables_come_before_the_usual_locations(bare_path, monkeypatch):
    usual = _fake_program(bare_path / "miniforge3" / "bin" / "mamba")
    from_var = _fake_program(bare_path / "elsewhere" / "conda")
    monkeypatch.setattr(command, "ENV_MANAGER_ROOTS", (str(bare_path / "miniforge3"),))
    monkeypatch.setenv("CONDA_EXE", from_var)
    assert command.locate_env_manager("mamba") == from_var
    mamba_exe = _fake_program(bare_path / "micro" / "mamba")
    monkeypatch.setenv("MAMBA_EXE", mamba_exe)
    assert command.locate_env_manager("mamba") == mamba_exe
    assert usual != mamba_exe


def test_the_configured_path_comes_first(bare_path, monkeypatch, tmp_path):
    configured = _fake_program(bare_path / "custom" / "mamba")
    monkeypatch.setenv("MAMBA_EXE", _fake_program(bare_path / "other" / "mamba"))
    config = tmp_path / "metawrap2.toml"
    config.write_text('[settings]\nenv_manager = "%s"\n' % configured)
    settings = load_settings(str(config))
    assert settings.env_manager == configured and "env_manager" in SETTINGS_KEYS
    command.set_configured_env_manager(settings.env_manager)
    assert command.locate_env_manager("mamba") == configured


def test_a_manager_on_path_is_used_by_name(bare_path, monkeypatch):
    _fake_program(bare_path / "empty_bin" / "mamba")
    monkeypatch.setenv("MAMBA_EXE", _fake_program(bare_path / "other" / "mamba"))
    assert command.locate_env_manager("mamba") == "mamba"
    assert command.ensure_env_manager() and command.runner.env_manager == "mamba"


def test_a_file_that_is_not_executable_is_skipped(bare_path, monkeypatch):
    path = bare_path / "miniforge3" / "bin" / "mamba"
    path.parent.mkdir(parents=True)
    path.write_text("")
    os.chmod(path, 0o644)
    monkeypatch.setattr(command, "ENV_MANAGER_ROOTS", (str(bare_path / "miniforge3"),))
    assert command.locate_env_manager("mamba") is None
