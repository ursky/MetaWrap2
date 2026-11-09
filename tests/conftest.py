import os
import sys

import pytest

# Allow tests to import the package from src/ without an install step.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))


@pytest.fixture(autouse=True)
def _conda_env_policy(request, tmp_path_factory, monkeypatch):
    """Pin ``use_conda_envs`` for every test instead of inheriting the machine's config.

    Most tests fake the tool call and want no env wrapping. Tests about env resolution mark
    themselves ``@pytest.mark.conda_envs``. Without this, results would depend on the local
    ``~/.metawrap2/config.toml`` or ``$METAWRAP2_CONFIG``.
    """
    wanted = request.node.get_closest_marker("conda_envs") is not None
    path = tmp_path_factory.mktemp("config") / "metawrap2.toml"
    # The shared score cache is off, so no test reads or writes scores kept between runs. A test
    # that writes its own config without the key gets a temporary default directory instead.
    path.write_text(
        '[settings]\nuse_conda_envs = %s\nscore_cache_dir = ""\n' % ("true" if wanted else "false")
    )
    monkeypatch.setenv("METAWRAP2_CONFIG", str(path))
    from metawrap2 import scorecache

    monkeypatch.setattr(
        scorecache, "DEFAULT_DIR", str(tmp_path_factory.mktemp("score_cache_default"))
    )


#: Fake machine seen by every test: 32 cores, 128 GB, no GPU, 1.8 TB free of 2 TB. It is larger
#: than any limit a test exercises, so thread caps and disk ceilings are tested against the cap.
TEST_MACHINE_CORES = 32
TEST_MACHINE_MEMORY = 128 * 1024**3
TEST_DISK = (2 * 1024**4, 200 * 1024**3, 2 * 1024**4 - 200 * 1024**3)  # total, used, free


@pytest.fixture(autouse=True)
def _no_env_manager_off_path(monkeypatch):
    """Find mamba or conda on PATH only, as on a clean CI runner.

    Off PATH, the command runner also looks in $MAMBA_EXE, $CONDA_EXE and the usual install
    locations, which differ between machines. Tests of that search set them up themselves.
    """
    from metawrap2 import command

    monkeypatch.delenv("MAMBA_EXE", raising=False)
    monkeypatch.delenv("CONDA_EXE", raising=False)
    monkeypatch.setattr(command, "ENV_MANAGER_ROOTS", ())
    monkeypatch.setattr(command, "_configured_env_manager", "")


@pytest.fixture(autouse=True)
def _fixed_machine(request, monkeypatch):
    """Give every test the same cores, memory, GPUs and free disk, whatever the host.

    Several modules clamp to the machine (threads to cores, downloads to free disk), so without
    this the results would differ between a workstation and a CI runner.

    Tests that measure the real host mark themselves ``@pytest.mark.real_machine``. Tests that need
    a different machine monkeypatch on top of this, which takes precedence.
    """
    if request.node.get_closest_marker("real_machine") is not None:
        return
    import collections
    import shutil

    from metawrap2 import limits

    monkeypatch.setattr(limits, "machine_cores", lambda: TEST_MACHINE_CORES)
    monkeypatch.setattr(limits, "machine_memory_bytes", lambda: TEST_MACHINE_MEMORY)
    monkeypatch.setattr(limits, "machine_gpus", lambda: 0)
    usage = collections.namedtuple("usage", "total used free")(*TEST_DISK)
    monkeypatch.setattr(shutil, "disk_usage", lambda _path: usage)
