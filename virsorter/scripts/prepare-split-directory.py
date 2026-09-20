#!/usr/bin/env python3
"""Reset a checkpoint-owned split directory without deleting its previous contents.

Snakemake directory() previously removed this directory before a rerun. A
completion-file checkpoint must explicitly reset stale split membership instead.
Keep the old directory as a sibling recovery artifact; never follow a symlink.
"""
import os
from pathlib import Path
import sys
import uuid


def prepare(path, root=None):
    root = Path(root or os.getcwd()).resolve()
    target = Path(path)
    if not target.is_absolute():
        target = root / target
    if target.is_symlink() or root not in target.resolve().parents:
        raise ValueError('Split directory must be a real child of the working directory')
    if target.exists():
        if not target.is_dir():
            raise ValueError('Split directory path is not a directory')
        previous = target.with_name(target.name + '.previous-' + uuid.uuid4().hex)
        target.rename(previous)
    target.mkdir(parents=True, exist_ok=False)


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: prepare-split-directory.py SPLIT_DIRECTORY')
    prepare(sys.argv[1])
