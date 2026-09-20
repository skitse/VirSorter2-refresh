"""Snakemake 9 command construction; scientific tools stay in rule-level envs."""
import shlex
import base64
import json
from pathlib import Path
import subprocess


def incomplete_targets(directory):
    """Read Snakemake 9 file-persistence records, without modifying journals.

    Dynamic checkpoint producers may only enter the DAG after its initial
    incomplete check. Explicit forcerun files make the native DAG revisit them.
    Unknown records fail closed; this is not a second dependency graph.
    """
    root = Path(directory).absolute()
    if (root / '.snakemake' / 'metadata.db').exists():
        raise ValueError('Database-backed persistence is not supported by this recovery fork; use a new file-backed run')
    folder = root / '.snakemake' / 'incomplete'
    if not folder.exists():
        return []
    result = []
    for path in sorted(folder.rglob('*')):
        if path.is_symlink():
            raise ValueError('Symlink in native incomplete metadata; inspect run before resuming')
        if not path.is_file():
            continue
        parts = path.relative_to(folder).parts
        if any(not part.startswith('@') for part in parts[:-1]):
            raise ValueError('Unknown Snakemake incomplete metadata layout')
        encoded = ''.join(part[1:] for part in parts[:-1]) + parts[-1]
        try:
            target = base64.b64decode(encoded, altchars=b'-_', validate=True).decode('utf-8')
            record = json.loads(path.read_text())
        except (ValueError, UnicodeError) as exc:
            raise ValueError('Invalid native incomplete record: ' + str(path)) from exc
        if not isinstance(record, dict) or 'external_jobid' not in record:
            raise ValueError('Unknown native incomplete record schema')
        if record['external_jobid'] is not None:
            raise ValueError('External job owns incomplete output; verify it is stopped before recovery: ' + target)
        destination = Path(target)
        resolved = (destination if destination.is_absolute() else root / destination).resolve()
        if not resolved.is_relative_to(root.resolve()) or resolved == root.resolve():
            raise ValueError('Incomplete output is outside the requested workdir')
        if not target or '\x00' in target or target.startswith('-'):
            raise ValueError('Invalid incomplete output identity')
        result.append(target)
    return result


def command(snakefile, directory, jobs, *, configfile=None, config=None,
            conda_prefix=None, use_conda=True, profile=None, dryrun=False,
            verbose=False, targets=(), force=(), extra=()):
    extra = list(extra)
    for argument in extra:
        if argument.startswith('--forcer') or argument.startswith('-R'):
            raise ValueError('Forwarded --forcerun/-R can override incomplete repair; finish recovery before using custom force targets')
        if argument.startswith('--persistence-'):
            raise ValueError('This recovery fork requires native file persistence; profile backend settings are overridden explicitly')
    if int(jobs) < 1:
        raise ValueError('jobs must be positive')
    argv = ['snakemake', '--snakefile', str(snakefile), '--directory', str(directory),
            '--jobs', str(jobs), '--rerun-incomplete', '--latency-wait', '600',
            '--persistence-backend', 'file']
    if profile:
        argv += ['--profile', str(profile)]
    else:
        argv += ['--cores', str(jobs)]
    if configfile:
        argv += ['--configfile', str(configfile)]
    if config:
        argv += ['--config'] + [f'{key}={value}' for key, value in config.items()]
    if use_conda:
        argv += ['--software-deployment-method', 'conda', '--conda-frontend', 'conda']
        if conda_prefix:
            argv += ['--conda-prefix', str(conda_prefix)]
    if dryrun:
        argv.append('--dry-run')
    if verbose:
        argv.append('--printshellcmds')
    repair = [] if '--unlock' in extra else incomplete_targets(directory)
    forced = list(dict.fromkeys([*map(str, force), *repair]))
    if forced:
        argv += ['--forcerun', *forced]
    # Preserve caller token boundaries. No shell parsing, interpolation or globbing.
    extra = list(extra)
    if extra[:1] == ['--']:
        extra = extra[1:]
    argv += extra
    if targets:
        argv += ['--', *map(str, targets)]
    return argv


def execute(argv):
    import logging
    logging.info('Executing: %s', shlex.join(argv))
    return subprocess.run(argv, check=True)
