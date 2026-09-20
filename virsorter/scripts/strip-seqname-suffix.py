#!/usr/bin/env python3
"""Portable, byte-preserving equivalents of finalize's suffix-off sed commands.

Modes retain the original boundaries and first-match-per-line behavior, not a
FASTA parser's normalization. In particular, end-of-line means before LF only:
a trailing CR or description prevents an end-anchored match, just as in sed.
Score whitespace is the ASCII/POSIX C-locale whitespace used by these tables.
"""

import argparse
import os
import re
import stat
import tempfile


_SUFFIX = rb'(?:full|[0-9]+_partial|lt2gene)'
_TRANSFORMS = {
    'score': ((re.compile(rb'\|\|' + _SUFFIX + rb'([ \t\r\v\f]+)'), rb'\1'),),
    'sequence': ((re.compile(rb'\|\|' + _SUFFIX + rb'\Z'), b''),),
    'dramv-affi': (
        (re.compile(rb'__' + _SUFFIX + rb'(\|[0-9]+\|[cl]\Z)'), rb'\1'),
        (re.compile(rb'__' + _SUFFIX + rb'(__[0-9]+\|)'), rb'\1'),
    ),
    'dramv-fasta': ((re.compile(rb'__' + _SUFFIX + rb'(-cat_[1-6]\Z)'), rb'\1'),),
}


def strip_line(line, mode):
    """Transform one binary line, preserving its optional LF byte."""
    ending = b'\n' if line.endswith(b'\n') else b''
    body = line[:-1] if ending else line
    for pattern, replacement in _TRANSFORMS[mode]:
        body = pattern.sub(replacement, body, count=1)
    return body + ending


def strip_file(path, mode):
    """Stream into a sibling temporary file and atomically replace one file.

    Preserve permission bits. Failures before replacement leave the input
    untouched and remove the temporary file. Multiple input files are replaced
    independently, not as a transaction. Like sed -i, replace a symlink itself.
    """
    path = os.path.abspath(path)
    temporary = None
    try:
        with open(path, 'rb') as source:
            permissions = stat.S_IMODE(os.fstat(source.fileno()).st_mode)
            with tempfile.NamedTemporaryFile(
                mode='wb', dir=os.path.dirname(path),
                prefix='.' + os.path.basename(path) + '.', delete=False,
            ) as target:
                temporary = target.name
                for line in source:
                    target.write(strip_line(line, mode))
                target.flush()
                os.fchmod(target.fileno(), permissions)
                os.fsync(target.fileno())
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            os.unlink(temporary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode', choices=sorted(_TRANSFORMS))
    parser.add_argument('files', nargs='+')
    args = parser.parse_args()
    for path in args.files:
        strip_file(path, args.mode)


if __name__ == '__main__':
    main()
