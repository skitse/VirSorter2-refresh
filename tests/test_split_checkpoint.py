import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('prepare_split',ROOT/'virsorter/scripts/prepare-split-directory.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


def test_rebuild_does_not_leave_stale_splits_or_delete_evidence(tmp_path):
    old=tmp_path/'iter-0/split';old.mkdir(parents=True)
    (old/'old-extra.part').write_text('prior evidence')
    m.prepare('iter-0/split',tmp_path)
    assert list(old.iterdir())==[]
    backups=list(old.parent.glob('split.previous-*'))
    assert len(backups)==1
    assert (backups[0]/'old-extra.part').read_text()=='prior evidence'


def test_cannot_reset_root_external_or_symlink(tmp_path):
    outside=tmp_path.parent/'outside-split-test'
    for candidate in (tmp_path,outside):
        with pytest.raises(ValueError):m.prepare(candidate,tmp_path)
    target=tmp_path/'real';target.mkdir();link=tmp_path/'link';link.symlink_to(target,target_is_directory=True)
    with pytest.raises(ValueError):m.prepare(link,tmp_path)
    assert target.is_dir()


def test_no_directory_checkpoint_owns_other_rule_outputs():
    for name in ('preprocess.smk','extract-feature.smk','classify.smk','train-feature.smk'):
        text=(ROOT/'virsorter/rules'/name).read_text()
        assert 'output: directory(' not in text
        assert 'prepare-split-directory.py' in text
        assert '/.split.done' in text
