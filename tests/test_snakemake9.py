"""Controller/CLI smoke and synthetic resume tests, no biological databases."""
import os
from pathlib import Path
import signal
import subprocess
import sys
import time

from click.testing import CliRunner
import pytest

from virsorter.snakemake_runner import command
from virsorter.virsorter import cli

ROOT=Path(__file__).resolve().parents[1]


def test_command_keeps_lock_and_token_boundaries():
    argv=command('/repo with space/Snakefile','/output with space',4,
                 configfile='/output with space/config.yaml',conda_prefix='/db/envs',
                 targets=['all'],extra=['--latency-wait','7'])
    assert '--nolock' not in argv
    assert '--conda-frontend' in argv and 'mamba' not in argv
    assert argv[argv.index('--snakefile')+1]=='/repo with space/Snakefile'
    assert argv[-2:]==['--','all']
    assert '--software-deployment-method' in argv


def test_bad_jobs_fail_before_launch():
    with pytest.raises(ValueError):command('Snakefile','out',0)


def test_external_profile_does_not_force_local_cores():
    argv=command('Snakefile','out',12,profile='cluster',use_conda=False)
    assert '--cores' not in argv
    assert '--software-deployment-method' not in argv
    assert '--profile' in argv


@pytest.mark.parametrize('provirus_off',[False,True])
def test_real_cli_dag_and_dryrun_preserves_metadata(tmp_path,provirus_off):
    db=tmp_path/'db'
    for group in ('dsDNAphage','ssDNA'):(db/'group'/group).mkdir(parents=True)
    out=tmp_path/'out'
    sentinel=out/'.snakemake'/'keep'
    sentinel.parent.mkdir(parents=True);sentinel.write_text('provenance')
    argv=['run','-w',str(out),'-i',str(ROOT/'test/8seq.fa'),'--db-dir',str(db),
          '--include-groups','dsDNAphage,ssDNA','--use-conda-off','--dryrun','-j','2']
    if provirus_off:argv+=['--provirus-off']
    argv+=['all']
    result=CliRunner().invoke(cli,argv)
    assert result.exit_code==0,repr(result.exception)+'\n'+result.output
    assert sentinel.read_text()=='provenance'


SNAKE=r'''
import os
wildcard_constraints:
    g='a|b',
    i='[0-9]+'
rule all:
    input: expand('merged/{g}.txt', g=['a','b'])
checkpoint split:
    output: touch('split/{g}/.split.done')
    params: directory='split/{g}'
    shell: 'mkdir -p {params.directory}; echo source > {params.directory}/0.part'
rule produce:
    input: 'split/{g}/{i}.part'
    output: temp('split/{g}/{i}.part.result')
    shell: 'printf partial > {output}; touch started-{wildcards.g}; sleep 2; echo completed > {output}'
def sources(wc):
    folder=os.path.dirname(checkpoints.split.get(g=wc.g).output[0])
    ids=glob_wildcards(os.path.join(folder,'{i}.part')).i
    return expand(os.path.join(folder,'{i}.part.result'),i=ids)
localrules: merge
rule merge:
    input: sources
    output: 'merged/{g}.txt'
    shell: 'mkdir -p merged; cat {input} > {output}'
'''


@pytest.mark.parametrize('stop_signal',[signal.SIGTERM,signal.SIGKILL])
def test_real_interrupted_checkpoint_resume(tmp_path,stop_signal):
    snake=tmp_path/'Snakefile';snake.write_text(SNAKE)
    args=['snakemake','-s',str(snake),'--directory',str(tmp_path),'--cores','4','--notemp','--rerun-incomplete']
    first=subprocess.run(args+['all'],capture_output=True,text=True,timeout=45)
    assert first.returncode==0,first.stdout+first.stderr
    for p in tmp_path.glob('started-*'):p.unlink()
    interrupted_log=tmp_path/'interrupted.log'
    with interrupted_log.open('w') as handle:
        proc=subprocess.Popen(args+['--forcerun','produce','--','all'],stdout=handle,stderr=subprocess.STDOUT,start_new_session=True)
        try:
            deadline=time.monotonic()+30
            while time.monotonic()<deadline and not list(tmp_path.glob('started-*')):
                assert proc.poll() is None,interrupted_log.read_text()
                time.sleep(.02)
            assert list(tmp_path.glob('started-*')),interrupted_log.read_text()
            os.killpg(proc.pid,stop_signal)
            proc.wait(timeout=20)
        finally:
            if proc.poll() is None:
                os.killpg(proc.pid,signal.SIGKILL);proc.wait()
    if stop_signal == signal.SIGKILL:
        # Exact regression: producer output still exists but is incomplete;
        # the old engine could schedule merge while replacing this partial file.
        assert any(p.read_text() == 'partial' for p in tmp_path.glob('split/*/*.result'))
    # Only this test owns this directory; the killed process group has exited.
    unlock=subprocess.run(args+['--unlock'],capture_output=True,text=True,timeout=30)
    assert unlock.returncode==0,unlock.stdout+unlock.stderr
    resume_args=command(snake,tmp_path,4,use_conda=False,targets=['all'],extra=['--notemp'])
    if stop_signal == signal.SIGKILL:
        assert '--forcerun' in resume_args
    resumed=subprocess.run(resume_args,capture_output=True,text=True,timeout=45)
    assert resumed.returncode==0,resumed.stdout+resumed.stderr
    for group in ('a','b'):
        assert (tmp_path/'merged'/f'{group}.txt').read_text()=='completed\n'
        producer=tmp_path/'split'/group/'0.part.result'
        # Graceful SIGTERM can remove temp outputs while retaining the prior
        # valid final result. A surviving incomplete output must be rebuilt.
        if stop_signal == signal.SIGKILL or producer.exists():
            assert producer.read_text()=='completed\n'
    assert 'No such file' not in resumed.stdout+resumed.stderr


def test_native_incomplete_record_becomes_explicit_force_target(tmp_path):
    import base64,json
    from virsorter.snakemake_runner import incomplete_targets
    folder=tmp_path/'.snakemake/incomplete';folder.mkdir(parents=True)
    target='split/a/result.tsv'
    record=folder/base64.urlsafe_b64encode(target.encode()).decode()
    record.write_text(json.dumps({'external_jobid':None}))
    assert incomplete_targets(tmp_path)==[target]
    argv=command('Snakefile',tmp_path,2,use_conda=False,targets=['all'])
    assert argv[argv.index('--forcerun')+1]==target
    record.write_text(json.dumps({'external_jobid':'active-cluster-job'}))
    with pytest.raises(ValueError,match='External job'):
        incomplete_targets(tmp_path)


def test_native_incomplete_external_path_is_rejected(tmp_path):
    import base64,json
    from virsorter.snakemake_runner import incomplete_targets
    folder=tmp_path/'.snakemake/incomplete';folder.mkdir(parents=True)
    (folder/base64.urlsafe_b64encode(b'../foreign/result').decode()).write_text(json.dumps({'external_jobid':None}))
    with pytest.raises(ValueError,match='outside'):
        incomplete_targets(tmp_path)


def test_legacy_run_requires_new_workdir_without_mutation(tmp_path):
    out=tmp_path/'legacy';out.mkdir()
    old='legacy scientific settings\n'
    (out/'config.yaml').write_text(old)
    result=CliRunner().invoke(cli,['run','-w',str(out),'-i',str(ROOT/'test/8seq.fa'),'--dryrun'])
    assert result.exit_code!=0
    assert 'in-place upgrade' in result.output
    assert (out/'config.yaml').read_text()==old
    assert not (out/'.virsorter-layout.json').exists()


@pytest.mark.parametrize('override',[['--forcerun','other_rule'],['-R','other_rule'],['--forcerun=other_rule'],['-Rother_rule']])
def test_forwarded_force_cannot_erase_repair_targets(tmp_path,override):
    with pytest.raises(ValueError,match='override incomplete repair'):
        command('Snakefile',tmp_path,2,use_conda=False,extra=override)


def test_native_parser_receives_file_backend_and_all_repair_targets(tmp_path):
    import base64,json
    from snakemake.cli import get_argument_parser
    folder=tmp_path/'.snakemake/incomplete';folder.mkdir(parents=True)
    target='split/a/result.tsv'
    (folder/base64.urlsafe_b64encode(target.encode()).decode()).write_text(json.dumps({'external_jobid':None}))
    argv=command('Snakefile',tmp_path,2,profile='cluster',use_conda=False,force=['classify'],targets=['all'])
    parsed=get_argument_parser().parse_args(argv[1:])
    assert parsed.persistence_backend.name=='FILE'
    assert set(parsed.forcerun)=={'classify',target}
    with pytest.raises(ValueError,match='file persistence'):
        command('Snakefile',tmp_path,2,extra=['--persistence-backend','db'])


def test_existing_database_persistence_not_silently_ignored(tmp_path):
    from virsorter.snakemake_runner import incomplete_targets
    (tmp_path/'.snakemake').mkdir()
    (tmp_path/'.snakemake/metadata.db').touch()
    with pytest.raises(ValueError,match='Database-backed'):
        incomplete_targets(tmp_path)
