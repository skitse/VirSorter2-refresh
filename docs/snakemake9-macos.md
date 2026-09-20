# Snakemake 9 and macOS compatibility branch

This is an experimental fork of VirSorter2 2.2.4. It is not an upstream release,
and it does not claim validated end-to-end prediction equivalence yet.

## Scope

- Python 3.11+ controller with Snakemake >=9.23.1,<10.
- Updated configuration import and workflow source-directory API.
- Argument-vector CLI execution, native directory locking, retained `.snakemake`
  recovery metadata, and non-destructive dry-run cleanup behavior.
- Named checkpoint inputs retained even before checkpoints resolve in dry-runs.
- Split checkpoints own completion files, not directories containing outputs of
  other rules. Rebuilding retains previous split directories as recovery artifacts.
- Native Snakemake 9 file-backend incomplete journals are validated and their
  outputs supplied as explicit native `--forcerun` targets. Dynamic checkpoint
  expansion otherwise can miss existing incomplete producers even on 9.23.1.
  Unknown journal formats, external-job ownership and paths outside the workdir
  are refused rather than guessed. File persistence is explicit, overriding
  profile defaults; database-backed runs need a new file-backed directory.
  Forwarded `--forcerun` / `-R` and persistence overrides are rejected so they
  cannot silently replace mandatory recovery targets.
- Portable, atomic Python suffix transformations replace GNU-only `sed -i -E`.
- Existing scientific thresholds, trained models, and rule-level environment
  (including the legacy scikit-learn version) are unchanged.

The old controller can schedule a consumer prematurely when a completed
checkpoint is reused but its producer must be rerun as incomplete. Tests cover
an actual interrupted producer process group and a subsequent native resume.
Do not fix this by deleting metadata or ignoring missing results.

## Installation

```bash
python3.12 -m venv .venv
. .venv/bin/activate
python -m pip install -e . pytest
virsorter --help
python -m pytest -q tests
```

For production prediction, the controller also needs a working Conda installation
and the original databases and rule-level environments. `vs2-external-deps.yaml`
now describes the controller only. Do not install modern scikit-learn over the
legacy model environment merely to satisfy the controller.

## macOS support boundary

| Component | Evidence / status |
|---|---|
| CLI help, config generation, main DAG dry-run | Locally tested on macOS arm64 |
| Snakemake 9 synthetic checkpoint interruption/resume | Locally tested |
| Suffix transformations vs original sed expressions | Locally tested with BSD sed; Linux CI also exercises sed |
| Linux/macOS Python 3.11/3.12 controller matrix | Defined in CI; inspect actual run status |
| Full prediction / database setup on Linux | Opt-in manual CI; not established by controller tests |
| Scientific dependencies and old serialized models on Apple Silicon | **Not validated** |
| Native macOS training or full model equivalence | **Not validated** |

Old scientific dependencies may not have native Apple Silicon packages. Do not
interpret a successful controller test as proof that those dependencies solve or
that model predictions match Linux. For now use the validated scientific platform
for real analyses; this branch remains isolated from existing production runs.

## Safe resume

Start this fork in a **new working directory**. It uses a different checkpoint
completion layout and refuses in-place adoption of an old run without its layout
receipt. Do not point it at an active production run.

After creating a run with this fork, resume with the same input, database,
scientific settings and working directory.
The CLI passes `--rerun-incomplete` and preserves native Snakemake locks. Never run
two controllers against one directory. Only unlock after proving the prior writer
and its children are gone. Input/database/scientific-setting changes need a new
work directory or an explicitly planned invalidation—not blind reuse.

The tests unlock **only their own temporary directory after terminating and
waiting for their test process group**. That is not a production auto-unlock policy.

## Known limitations

- CLI argument handling is token-safe, but some original rule shell snippets still
  assume shell-safe paths. Avoid spaces and shell metacharacters in production
  input, database and output paths until those snippets receive a full quoting audit.
- Dry-runs still create/update the native `config.yaml` and fork layout receipt; they do not delete
  `.snakemake` or temporary scientific results.
- A full fresh-run versus interrupted-run native-output comparison remains a
  release gate. No scientific equivalence claim is made from synthetic tests.
- Database-download/setup and custom training need separate real-tool acceptance.
