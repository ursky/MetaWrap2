"""Bin scores kept between runs, so an identical bin is not scored twice.

:class:`metawrap2.checkm.BinScoreCache` keeps the CheckM scores of one run under its output
directory, for ``--resume``. This module adds a second store shared by all of a user's runs. A bin
with the same sequences, scored by the same scorer, scorer version, database and options, reuses
the stored score: for example ``bin_refinement`` re-run with other ``-c``/``-x`` thresholds, or the
final bins scored again with ``score``. On one refinement, rescoring an identical candidate pool
took 1,934 s without such a cache and 9 s with it.

Layout: one JSON file per scorer configuration under ``[settings] score_cache_dir`` (default
``~/.metawrap2/score_cache``), named ``<scorer>.<hash of the configuration>.json``. It holds the
configuration and, for each bin, the six scored columns of a ``.stats.tsv`` row (completeness
through size), keyed by a hash of the bin's sequences (:func:`content_identity`). The per-run cache
keys on contig names, which identify a contig within one assembly but not across assemblies.

The configuration (:func:`scorer_configuration`) is the scorer (``checkm2``, ``checkm1_reduced``
or ``checkm1``), the scorer's conda package version read from its environment's ``conda-meta``,
a fingerprint of its database (:func:`database_identity`), CheckM2's ``--ttable`` choice and
:data:`SCHEMA`. If any part cannot be determined the run does not use the shared cache.

Several runs can share the directory. Files are replaced whole (temp file, then ``os.replace``),
so a reader never sees a partial file and does not lock. A writer re-reads the file and merges its
new scores under an ``fcntl`` lock on ``<file>.lock``. If the lock is not free within
:data:`LOCK_TIMEOUT_SECONDS` the write is skipped, and the scores are kept for the next attempt.
A failure to read or write the cache never fails a run.

CheckM2's DIAMOND hits for a protein can depend on the other proteins in the same call, so a score
reused from another run can differ in the last digits from a fresh score of the same bin in a new
grouping. The per-run cache accepts the same difference between bin sets of one run.
"""

from __future__ import annotations

import errno
import fcntl
import glob
import hashlib
import json
import os
import shutil
import threading
import time
from typing import Any, Callable, Dict, List, Optional, Tuple

#: Version of the stored row format. Raise it when the ``.stats.tsv`` columns or their formatting
#: change, so scores written by older code are not reused.
SCHEMA = 1

#: Where the cache lives when ``[settings] score_cache_dir`` is not set: the user's MetaWrap2
#: directory, beside ``config.toml`` and ``history.jsonl``.
DEFAULT_DIR = os.path.join("~", ".metawrap2", "score_cache")

#: How long a run waits for another run's write before it skips its own.
LOCK_TIMEOUT_SECONDS = 10.0
_LOCK_POLL_SECONDS = 0.1

#: Columns stored per bin: STATS_HEADER[1:7], completeness to size. The bin name and an optional
#: binner column are filled in by the run that reuses the score.
N_FIELDS = 6

#: Scorer name -> (conda package, executable, database key, label for messages).
SCORERS: Dict[str, Tuple[str, str, str, str]] = {
    "checkm2": ("checkm2", "checkm2", "CHECKM2DB", "CheckM2"),
    "checkm1_reduced": ("checkm-genome", "checkm", "CHECKM_DB", "CheckM1"),
    "checkm1": ("checkm-genome", "checkm", "CHECKM_DB", "CheckM1"),
}

#: Key under which a run's use of the shared cache is written to ``run_config.json``.
RECORD_KEY = "score_cache"


def scorer_name(checkm2: bool, full_tree: bool = False) -> str:
    """``checkm2``, ``checkm1_reduced`` or ``checkm1`` (full tree), as in bin_refinement."""
    if checkm2:
        return "checkm2"
    return "checkm1" if full_tree else "checkm1_reduced"


def cache_dir(setting: Optional[str]) -> Optional[str]:
    """The shared cache directory for a ``score_cache_dir`` value, or None if it is turned off.

    None (key not set) means :data:`DEFAULT_DIR`; an empty string turns the cache off.
    """
    if setting is None:
        setting = DEFAULT_DIR
    setting = setting.strip()
    if not setting:
        return None
    return os.path.abspath(os.path.expanduser(setting))


def content_identity(path: str) -> str:
    """A key for a bin's contents: a hash of its sequences, independent of order and names."""
    from .io.seqio import iter_fasta

    digests = sorted(hashlib.sha1(seq.encode()).hexdigest() for _header, seq in iter_fasta(path))
    return hashlib.sha256("\n".join(digests).encode()).hexdigest()


def database_identity(key: str, path: str) -> Optional[str]:
    """A fingerprint of the database at *path* (a file or a directory), or None if absent.

    A file is identified by its name, size and mtime; a directory by
    :func:`metawrap2.dbcheck.fingerprint`, which does the same for its defining files.
    """
    if not path:
        return None
    if os.path.isfile(path):
        try:
            stat = os.stat(path)
        except OSError:
            return None
        text = "%s:%d:%d" % (os.path.basename(path), stat.st_size, int(stat.st_mtime))
        return hashlib.sha256(text.encode()).hexdigest()[:16]
    if os.path.isdir(path):
        from .dbcheck import fingerprint

        return fingerprint(path, key)
    return None


def _env_root(env: Optional[str], executable: str) -> Optional[str]:
    """The conda prefix the scorer runs from: the named env, or the one holding it on PATH."""
    if env:
        from .config import env_prefix

        return env_prefix(env)
    found = shutil.which(executable)
    if not found:
        return None
    return os.path.dirname(os.path.dirname(os.path.realpath(found)))


def tool_version(package: str, executable: str, env: Optional[str]) -> Optional[str]:
    """``package=version=build`` from the scorer env's ``conda-meta``, or None if not found.

    Read from the package record, so the tool is not run.
    """
    root = _env_root(env, executable)
    if not root:
        return None
    found = []
    for record in glob.glob(os.path.join(root, "conda-meta", package + "-*.json")):
        parts = os.path.basename(record)[: -len(".json")].rsplit("-", 2)
        if len(parts) == 3 and parts[0] == package:
            found.append("=".join(parts))
    if len(found) != 1:
        return None
    return found[0]


def scorer_configuration(
    databases: Dict[str, str], scorer: str, ttable: Optional[int], env: Optional[str]
) -> Tuple[Optional[Dict[str, Any]], str]:
    """``(configuration, "")``, or ``(None, reason)`` if part of it cannot be determined.

    *databases* is the config's ``[databases]`` table, *env* the conda env the scorer runs in
    (None when tools run from PATH).
    """
    package, executable, db_key, label = SCORERS[scorer]
    path = databases.get(db_key, "")
    if not path:
        return None, "%s is not set, so the database %s uses is not known" % (db_key, label)
    database = database_identity(db_key, path)
    if database is None:
        return None, "the %s database at %s could not be read" % (label, path)
    version = tool_version(package, executable, env)
    if version is None:
        return None, "the installed %s version could not be read from its conda environment" % (
            label
        )
    config: Dict[str, Any] = {
        "schema": SCHEMA,
        "scorer": scorer,
        "tool": version,
        "database": database,
    }
    if scorer == "checkm2":
        # CheckM2 tries tables 11 and 4 unless --ttable picks one.
        config["ttable"] = str(ttable) if ttable else "both"
    return config, ""


def _acquire(handle, timeout: float) -> bool:
    """Take an exclusive lock on *handle*, waiting at most *timeout* seconds."""
    deadline = time.monotonic() + timeout
    while True:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return True
        except OSError as exc:
            if exc.errno not in (errno.EAGAIN, errno.EACCES, errno.EWOULDBLOCK):
                raise
        if time.monotonic() >= deadline:
            return False
        time.sleep(_LOCK_POLL_SECONDS)


class SharedScoreCache:
    """The stored scores for one scorer configuration, and this run's use of them."""

    def __init__(
        self,
        directory: str,
        config: Dict[str, Any],
        say: Optional[Callable[[str], None]] = None,
        lock_timeout: float = LOCK_TIMEOUT_SECONDS,
    ) -> None:
        self.directory = directory
        self.config = dict(config)
        digest = hashlib.sha256(json.dumps(self.config, sort_keys=True).encode()).hexdigest()
        self.path = os.path.join(directory, "%s.%s.json" % (self.config["scorer"], digest[:16]))
        self.lock_path = self.path + ".lock"
        self.lock_timeout = lock_timeout
        self.say = say or (lambda _msg: None)
        self._lock = threading.RLock()
        self._flush_lock = threading.Lock()
        self._warned_write = False
        self.rows: Dict[str, List[str]] = self._read()
        self.pending: Dict[str, List[str]] = {}
        self.loaded = len(self.rows)
        self.hits = self.misses = self.stored = self.writes_skipped = 0

    def _read(self) -> Dict[str, List[str]]:
        """The stored rows, or {} if the file is missing, unreadable or for another configuration."""
        try:
            with open(self.path) as fh:
                data = json.load(fh)
        except (OSError, ValueError):
            return {}
        if not isinstance(data, dict) or data.get("config") != self.config:
            return {}
        rows = data.get("rows")
        if not isinstance(rows, dict):
            return {}
        return {
            str(key): [str(v) for v in value]
            for key, value in rows.items()
            if isinstance(value, list) and len(value) == N_FIELDS
        }

    def get(self, key: str) -> Optional[List[str]]:
        """The stored columns for content key *key*, counting a hit or a miss."""
        with self._lock:
            found = self.pending.get(key) or self.rows.get(key)
            if found is None:
                self.misses += 1
                return None
            self.hits += 1
            return list(found)

    def put(self, key: str, values: List[str]) -> None:
        """Queue a newly scored bin's columns for the next :meth:`flush`."""
        if len(values) != N_FIELDS:
            return
        with self._lock:
            if key not in self.rows:
                self.pending[key] = list(values)

    def flush(self) -> bool:
        """Merge queued scores into the file. False if they were not written (they stay queued)."""
        with self._flush_lock:
            with self._lock:
                pending = dict(self.pending)
            if not pending:
                return True
            try:
                os.makedirs(self.directory, exist_ok=True)
                with open(self.lock_path, "a") as handle:
                    if not _acquire(handle, self.lock_timeout):
                        self.writes_skipped += 1
                        self.say(
                            "The shared score cache is locked by another run; %d new bin score(s) "
                            "were not added to it this time." % len(pending)
                        )
                        return False
                    try:
                        current = self._read()
                        added = 0
                        for key, values in pending.items():
                            if key not in current:
                                current[key] = values
                                added += 1
                        payload = json.dumps(
                            {"config": self.config, "rows": current}, separators=(",", ":")
                        )
                        tmp = "%s.%d.%d.tmp" % (self.path, os.getpid(), threading.get_ident())
                        with open(tmp, "w") as fh:
                            fh.write(payload)
                        os.replace(tmp, self.path)
                    finally:
                        fcntl.flock(handle, fcntl.LOCK_UN)
            except OSError as exc:
                if not self._warned_write:
                    self._warned_write = True
                    self.say(
                        "Could not write the shared score cache %s (%s); this run's scores are "
                        "kept in its own output directory only." % (self.path, exc.strerror or exc)
                    )
                return False
            with self._lock:
                for key in pending:
                    self.pending.pop(key, None)
                self.rows.update(current)
                self.stored += added
            return True

    def summary(self) -> str:
        with self._lock:
            return "%d bin(s) reused from earlier runs, %d not found there, %d added" % (
                self.hits,
                self.misses,
                self.stored,
            )

    def as_record(self) -> Dict[str, Any]:
        """What ``run_config.json`` records about this run's use of the shared cache."""
        with self._lock:
            return {
                "status": "used",
                "file": self.path,
                "configuration": dict(self.config),
                "bins_stored_at_start": self.loaded,
                "hits": self.hits,
                "misses": self.misses,
                "added": self.stored,
                "writes_skipped_on_lock_timeout": self.writes_skipped,
                "not_written": len(self.pending),
            }


def _record(value: Any) -> None:
    """Attach *value* (or a callable giving it) to the active run's ``run_config.json``."""
    from .command import runner

    recorder = runner.recorder
    extra = getattr(recorder, "extra", None)
    if isinstance(extra, dict):
        extra[RECORD_KEY] = value


def open_for_run(
    setting: Optional[str],
    databases: Dict[str, str],
    *,
    enabled: bool,
    scorer: str,
    ttable: Optional[int] = None,
    env: Optional[str] = None,
    say: Optional[Callable[[str], None]] = None,
) -> Optional[SharedScoreCache]:
    """The shared cache for this run's scorer, or None (said once, and recorded) if not used.

    *setting* is ``[settings] score_cache_dir``; *enabled* is False under ``--no-score-cache``.
    A dry run never uses it.
    """
    from .command import runner

    log = say or (lambda _msg: None)
    if runner.dry_run:
        return None
    if not enabled:
        _record({"status": "off", "reason": "--no-score-cache"})
        log("Not using the shared score cache (--no-score-cache).")
        return None
    directory = cache_dir(setting)
    if directory is None:
        _record({"status": "off", "reason": "score_cache_dir is empty"})
        log("Not using the shared score cache (score_cache_dir is empty in the config).")
        return None
    config, reason = scorer_configuration(databases, scorer, ttable, env)
    if config is None:
        _record({"status": "off", "reason": reason})
        log("Not using the shared score cache: %s. Scores are kept for this run only." % reason)
        return None
    shared = SharedScoreCache(directory, config, say=log)
    _record(shared.as_record)
    log(
        "Shared score cache: %s (%d bin score(s) stored for this scorer, version and database)"
        % (shared.path, shared.loaded)
    )
    return shared
