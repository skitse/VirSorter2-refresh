"""Suffix-off parity tests; no databases, models, or workflow runtime needed."""

import importlib.util
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / 'virsorter' / 'scripts' / 'strip-seqname-suffix.py'
spec = importlib.util.spec_from_file_location('strip_suffix', SCRIPT)
helper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helper)

# Original expressions, copied verbatim from classify.smk (shell spelling).
SED = {
    'score': [r's/(\|\|full([[:space:]]+)|\|\|[0-9]+_partial([[:space:]]+)|\|\|lt2gene([[:space:]]+))/\2\3\4/;'],
    'sequence': [r's/(\|\|full$|\|\|[0-9]+_partial$|\|\|lt2gene$)//;'],
    'dramv-affi': [
        r's/(__full(\|[0-9]+\|(c|l)$)|__[0-9]+_partial(\|[0-9]+\|(c|l)$)|__lt2gene(\|[0-9]+\|(c|l)$))/\2\4\6/;',
        r's/(__full(__[0-9]+\|)|__[0-9]+_partial(__[0-9]+\|)|__lt2gene(__[0-9]+\|))/\2\3\4/;',
    ],
    'dramv-fasta': [r's/(__full(-cat_[1-6]$)|__[0-9]+_partial(-cat_[1-6]$)|__lt2gene(-cat_[1-6]$))/\2\3\4/;'],
}
SUFFIXES = (b'full', b'0_partial', b'12_partial', b'000_partial', b'lt2gene')


def examples(mode):
    """Explicit golden input/output pairs, including unusual legacy behavior."""
    yield b'', b''
    yield b'\xffunchanged\t|c\r\nACGT\x00\n', b'\xffunchanged\t|c\r\nACGT\x00\n'
    for suffix in SUFFIXES:
        for end in (b'\n', b''):
            if mode == 'score':
                for separator in (b'\t', b' ', b'\t  \t', b'\v', b'\f', b'\r'):
                    yield b'id||' + suffix + separator + b'0.9\t7' + end, b'id' + separator + b'0.9\t7' + end
                yield b'id||' + suffix + end, b'id||' + suffix + end
                yield b'id||' + suffix + b'\tother||full\t1' + end, b'id\tother||full\t1' + end
            elif mode == 'sequence':
                for prefix in (b'>circular|id', b'original\t1\t99\tid', b'\xffid'):
                    yield prefix + b'||' + suffix + end, prefix + end
                yield b'>id||' + suffix + b' description' + end, b'>id||' + suffix + b' description' + end
                yield b'>id||' + suffix + b'\r' + end, b'>id||' + suffix + b'\r' + end
                yield b'id||' + suffix + b'\t1' + end, b'id||' + suffix + b'\t1' + end
            elif mode == 'dramv-affi':
                for circular in (b'c', b'l'):
                    yield b'>id__' + suffix + b'|123|' + circular + end, b'>id|123|' + circular + end
                    yield b'gene__' + suffix + b'__01|x\t>id__' + suffix + b'|123|' + circular + end, b'gene__01|x\t>id|123|' + circular + end
                yield b'id__' + suffix + b'__42|\t0.9' + end, b'id__42|\t0.9' + end
                yield b'id__' + suffix + b'__1|id__full__2|' + end, b'id__1|id__full__2|' + end
                for tail in (b'|12|x', b'|12|c\t', b'|12|c\r', b'|x|l', b'__x|', b'__2\t', b'|2|c more'):
                    value = b'id__' + suffix + tail + end
                    yield value, value
            else:
                for category in range(1, 7):
                    tail = b'-cat_' + str(category).encode()
                    yield b'>id__' + suffix + tail + end, b'>id' + tail + end
                for tail in (b'-cat_0', b'-cat_7', b'-cat_10', b'-cat_1\r', b'-cat_1 desc', b'-cat_1\t'):
                    value = b'>id__' + suffix + tail + end
                    yield value, value
    for invalid in (b'partial', b'-1_partial', b'x_partial', b'FULL', b'lt2genes', b'fuller'):
        for tail in (b'', b'\t0.9', b'|1|c', b'__1|', b'-cat_1'):
            value = b'>id||' + invalid + tail + b'\n>id__' + invalid + tail + b'\n'
            yield value, value


class PortableSuffixTests(unittest.TestCase):
    def test_golden_bytes(self):
        for mode in SED:
            for source, expected in examples(mode):
                with self.subTest(mode=mode, source=source):
                    # Split on LF only, never bytes.splitlines (which treats CR).
                    with tempfile.TemporaryDirectory() as directory:
                        path = Path(directory) / 'data'
                        path.write_bytes(source)
                        helper.strip_file(path, mode)
                        self.assertEqual(path.read_bytes(), expected)

    def test_original_sed_parity(self):
        # BSD sed can validate expressions on macOS without using GNU -i.
        # GNU sed is also exercised if available (including Linux CI's sed).
        executables = list(dict.fromkeys(filter(None, (shutil.which('sed'), shutil.which('gsed')))))
        if not executables:
            self.skipTest('No sed available for independent expression parity')
        for executable in executables:
            for mode, expressions in SED.items():
                # Always LF-terminate oracle input: BSD sed adds a missing final
                # LF, unlike GNU sed -i; the golden test checks preservation.
                source = b''.join(value + (b'' if value.endswith(b'\n') else b'\n') for value, _ in examples(mode))
                command = [executable, '-E']
                for expression in expressions:
                    command.extend(['-e', expression])
                expected = subprocess.run(command, input=source, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, env=dict(os.environ, LC_ALL='C')).stdout
                actual = b''.join(helper.strip_line(line + b'\n', mode) for line in source.split(b'\n')[:-1])
                with self.subTest(sed=executable, mode=mode):
                    self.assertEqual(actual, expected)

    def test_cli_multiple_files_and_permissions(self):
        with tempfile.TemporaryDirectory() as directory:
            files = [Path(directory) / name for name in ('with space.fa', 'boundary.tsv')]
            for path in files:
                path.write_bytes(b'id||full\n\xffACGT\r\nlast||12_partial')
                path.chmod(0o640)
            subprocess.run([sys.executable, str(SCRIPT), 'sequence', *map(str, files)], check=True)
            for path in files:
                self.assertEqual(path.read_bytes(), b'id\n\xffACGT\r\nlast')
                self.assertEqual(stat.S_IMODE(path.stat().st_mode), 0o640)
            self.assertEqual(sorted(p.name for p in Path(directory).iterdir()), sorted(p.name for p in files))

    def test_replace_is_atomic_and_uses_sibling_temporary(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'data'
            original = b'>id||full\nACGT\n'
            expected = b'>id\nACGT\n'
            path.write_bytes(original)
            replace = os.replace

            def checked_replace(source, destination):
                self.assertEqual(Path(source).parent, path.parent)
                self.assertEqual(Path(destination), path)
                self.assertEqual(path.read_bytes(), original)
                self.assertEqual(Path(source).read_bytes(), expected)
                replace(source, destination)

            # An already-open reader must keep seeing the complete original.
            with path.open('rb') as old_reader:
                with mock.patch.object(helper.os, 'replace', side_effect=checked_replace) as call:
                    helper.strip_file(path, 'sequence')
                call.assert_called_once()
                self.assertEqual(old_reader.read(), original)
            self.assertEqual(path.read_bytes(), expected)

    def test_failure_keeps_original_and_removes_temporary(self):
        for target in ('strip_line', 'os.replace'):
            with self.subTest(target=target), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'data'
                original = b'>id||full\nACGT\n'
                path.write_bytes(original)
                parent, attribute = (helper.os, 'replace') if target == 'os.replace' else (helper, target)
                with mock.patch.object(parent, attribute, side_effect=OSError('injected failure')):
                    with self.assertRaises(OSError):
                        helper.strip_file(path, 'sequence')
                self.assertEqual(path.read_bytes(), original)
                self.assertEqual(list(Path(directory).iterdir()), [path])

    def test_both_finalize_branches_use_helper(self):
        text = (ROOT / 'virsorter' / 'rules' / 'classify.smk').read_text()
        self.assertNotIn('sed -i -E', text)
        self.assertEqual(text.count('strip-seqname-suffix.py score {output.score:q}'), 2)
        self.assertEqual(text.count('strip-seqname-suffix.py sequence {output.fa:q}'), 2)
        self.assertEqual(text.count('strip-seqname-suffix.py sequence {output.fa:q} {output.boundary:q}'), 1)
        self.assertEqual(text.count('strip-seqname-suffix.py dramv-affi '), 2)
        self.assertEqual(text.count('strip-seqname-suffix.py dramv-fasta '), 2)


if __name__ == '__main__':
    unittest.main()
