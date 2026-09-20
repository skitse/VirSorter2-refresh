"""Controller-side checks and explicit native zero-hit outcomes (no model changes)."""
import csv
import json
import math
from pathlib import Path
import re
import tempfile
import os
import uuid


def validate_request(workdir, workflow, min_score, min_length, jobs, label,
                     provirus_off, max_orf_per_seq, prep_for_dramv):
    if not math.isfinite(min_score) or not 0 <= min_score <= 1:
        raise ValueError('--min-score must be a finite number between 0 and 1')
    if min_length < 0 or jobs < 1:
        raise ValueError('--min-length must be nonnegative and --jobs must be positive')
    if label and not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', label):
        raise ValueError('--label must be a filename prefix, not a path or shell expression')
    if provirus_off and max_orf_per_seq != -1 and prep_for_dramv:
        raise ValueError('--max-orf-per-seq cannot be combined with --prep-for-dramv: subsampled ORFs are incompatible')
    if workflow == 'classify' and not (Path(workdir)/'config.yaml').is_file():
        raise ValueError('classify requires an existing completed feature-extraction run in the same -w directory; first run `virsorter run ... all`')


def preserve_previous_outcome(workdir, label=''):
    prefix = label + '-' if label else ''
    path = Path(workdir)/(prefix+'run-outcome.json')
    if path.exists():
        path.rename(path.with_name(path.name+'.previous-'+uuid.uuid4().hex))


def record_outcome(workdir, label=''):
    """Require native score/FASTA agreement; zero rows are valid, not failure."""
    root = Path(workdir)
    prefix = label + '-' if label else ''
    score = root/(prefix+'final-viral-score.tsv')
    fasta = root/(prefix+'final-viral-combined.fa')
    with score.open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'seqname', 'max_score'} <= set(reader.fieldnames or []):
            raise ValueError('Native score table is missing its required header')
        identifiers = set()
        for row in reader:
            name = row['seqname']
            if None in row or None in row.values() or not name or name in identifiers:
                raise ValueError('Malformed or duplicate native score record')
            identifiers.add(name)
    sequences = set()
    current = None
    length = 0
    with fasta.open() as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith('>'):
                if current is not None and length == 0:
                    raise ValueError('Empty native FASTA sequence')
                names = text[1:].split()
                if not names or names[0] in sequences:
                    raise ValueError('Malformed or duplicate native FASTA identifier')
                current = names[0]
                sequences.add(current)
                length = 0
            else:
                if current is None or set(text.upper()) - set('ACGTRYSWKMBDHVN'):
                    raise ValueError('Malformed native nucleotide FASTA')
                length += len(text)
    if current is not None and length == 0:
        raise ValueError('Empty native FASTA sequence')
    if identifiers != sequences:
        raise ValueError('Native score and FASTA sequence IDs disagree')
    result = {'status': 'completed_with_predictions' if identifiers else 'completed_no_viral_sequences',
              'predicted_sequences': len(identifiers),
              'scope': 'native score/FASTA output validation; not biological validation'}
    target = root/(prefix+'run-outcome.json')
    fd, temporary = tempfile.mkstemp(prefix='.outcome-', dir=root)
    try:
        with os.fdopen(fd, 'w') as handle:
            json.dump(result, handle, indent=2)
            handle.write('\n')
        os.replace(temporary, target)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)
    return result
