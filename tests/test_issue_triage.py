import json
from pathlib import Path
from click.testing import CliRunner
import pytest
from virsorter.virsorter import cli
from virsorter.run_validation import record_outcome

ROOT=Path(__file__).resolve().parents[1]

@pytest.mark.parametrize('options',[
    ['classify'],['--min-score','nan'],['--min-score','inf'],
    ['--label','../outside'],['--provirus-off','--max-orf-per-seq','20','--prep-for-dramv'],
])
def test_invalid_request_does_not_create_run(tmp_path,options):
    out=tmp_path/'new'
    result=CliRunner().invoke(cli,['run','-w',str(out),'-i',str(ROOT/'test/8seq.fa'),*options])
    assert result.exit_code!=0
    assert not out.exists()


def test_first_classify_has_actionable_message(tmp_path):
    result=CliRunner().invoke(cli,['run','-w',str(tmp_path/'new'),'-i',str(ROOT/'test/8seq.fa'),'classify'])
    assert 'first run' in result.output
    assert 'same -w' in result.output


def test_header_only_result_is_explicit_success(tmp_path):
    (tmp_path/'final-viral-score.tsv').write_text('seqname\tmax_score\n')
    (tmp_path/'final-viral-combined.fa').touch()
    assert record_outcome(tmp_path)['status']=='completed_no_viral_sequences'
    assert json.loads((tmp_path/'run-outcome.json').read_text())['predicted_sequences']==0


def test_nonzero_predictions_are_counted_with_label(tmp_path):
    (tmp_path/'run-final-viral-score.tsv').write_text('seqname\tmax_score\nx||lt2gene\tnan\n')
    (tmp_path/'run-final-viral-combined.fa').write_text('>x||lt2gene\nACGT\n')
    assert record_outcome(tmp_path,'run')['predicted_sequences']==1


@pytest.mark.parametrize('content',['','>different\nACGT\n','>x\n','ACGT\n'])
def test_missing_or_contradictory_results_do_not_get_success(tmp_path,content):
    (tmp_path/'final-viral-score.tsv').write_text('seqname\tmax_score\nx\t0.99\n')
    (tmp_path/'final-viral-combined.fa').write_text(content)
    with pytest.raises(ValueError):record_outcome(tmp_path)
    assert not (tmp_path/'run-outcome.json').exists()


def test_previous_success_is_retained_but_not_current_during_retry(tmp_path):
    from virsorter.run_validation import preserve_previous_outcome
    old=tmp_path/'run-outcome.json';old.write_text('{"status":"old success"}')
    preserve_previous_outcome(tmp_path)
    assert not old.exists()
    backups=list(tmp_path.glob('run-outcome.json.previous-*'))
    assert len(backups)==1
    assert json.loads(backups[0].read_text())['status']=='old success'
