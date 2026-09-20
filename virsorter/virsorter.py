import sys
import os
import logging
import multiprocessing
import subprocess
import glob
import shutil
import shlex
import json
import click

from snakemake.common.configfile import load_configfile
from virsorter.snakemake_runner import command as snakemake_command, execute as execute_snakemake
from ruamel.yaml import YAML
from virsorter import __version__
from virsorter.config import get_default_config, set_logger, make_config

set_logger()
logging.info(f'VirSorter {__version__}')
logging.info(' '.join(map(shlex.quote, sys.argv)))

def log_exception(msg):
    logging.critical(msg)
    logging.info('Documentation is available at: '
                    'https://github.com/jiarong/VirSorter2')
    logging.info('Issues can be raised at: '
                    'https://github.com/jiarong/VirSorter2/issues')
    logging.info('Or send an email to '
                    'virsorter2 at gmail.com if you do not use GitHub')
    sys.exit(1)

@click.group(context_settings=dict(help_option_names=["-h", "--help"]))
@click.version_option(__version__)
@click.pass_context
def cli(obj):
    """
    virsorter - workflow for identifying viral sequences
    """

def get_snakefile(f="Snakefile"):
    sf = os.path.join(os.path.dirname(os.path.abspath(__file__)), f)
    if not os.path.exists(sf):
        sys.exit("Unable to locate the Snakemake workflow file; tried %s" % sf)
    return sf

## run command

@cli.command(
    'run',
    context_settings=dict(ignore_unknown_options=True),
    short_help='run virsorter main workflow'
)
@click.argument(
    'workflow',
    default='all',
    type=click.Choice(['all', 'classify']),
#    show_default=True,
#    help='Execute only subworkflow.',
)
@click.option('-w',
    '--working-dir',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    help='output directory',
    default='.'
)
@click.option('-d',
    '--db-dir',
    required=False,
    type=click.Path(dir_okay=True,writable=False,resolve_path=True),
    help='database directory, default to the --db-dir set during installation',
)
@click.option('-i',
    '--seqfile',
    required=True,
    type=click.Path(resolve_path=True),
    help=('sequence file in fa or fq format (could be compressed '
            'by gzip or bz2)')
)
@click.option(
    '-l',
    '--label',
    default='',
    help=('add label as prefix to output files; this is useful when re-running '
            'classify with different filtering options'),
)
@click.option(
    '--include-groups',
    default='dsDNAphage,ssDNA',
    type=str,
    show_default=True,
    help=('classifiers of viral groups to included '
        '(comma separated and no space in between); available options are: '
        'dsDNAphage,NCLDV,RNA,ssDNA,lavidaviridae'),
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='max # of jobs allowed in parallel.',
)
@click.option(
    '--min-score',
    default=0.5,
    type=float,
    show_default=True,
    help='minimal score to be identified as viral',
)
@click.option(
    '--min-length',
    default=0,
    type=int,
    show_default=True,
    help=('minimal seq length required; all seqs shorter than this will '
            'be removed'),
)
@click.option(
    '--keep-original-seq',
    default=False,
    is_flag=True,
    show_default=True,
    help=('keep the original sequences instead of trimmed; By default, '
            'the untranslated regions at both ends of identified viral '
            'seqs are trimmed; circular sequences are modified to remove '
            'overlap between both ends and adjusted for the gene splitted '
            'into two ends;'),
)
@click.option(
    '--exclude-lt2gene',
    default=False,
    is_flag=True,
    show_default=True,
    help=('short seqs (less than 2 genes) does not have any scores, but '
            'those with hallmark genes are included as viral by default; ' 
            'use this option to exclude them'),
)
@click.option(
    '--prep-for-dramv',
    default=False,
    is_flag=True,
    show_default=True,
    help=('turn on generating viral seqfile and viral-affi-contigs.tab '
            'for DRAMv'),
)
@click.option(
    '--high-confidence-only',
    default=False,
    is_flag=True,
    show_default=True,
    help=('only output high confidence viral sequences; this is equivalent '
            'to screening final-viral-score.tsv with the following criteria: '
            '(max_score >= 0.9) OR (max_score >=0.7 AND hallmark >= 1)'),
)
@click.option(
    '--hallmark-required',
    default=False,
    is_flag=True,
    show_default=True,
    help='require hallmark gene on all viral seqs',
)
@click.option(
    '--hallmark-required-on-short',
    default=False,
    is_flag=True,
    show_default=True,
    help=('require hallmark gene on short viral seqs (length cutoff '
            'as "short" were set by "MIN_SIZE_ALLOWED_WO_HALLMARK_GENE" in '
            'template-config.yaml file, default 3kbp); this can reduce '
            'false positives at reasonable cost of sensitivity'),
)
@click.option(
    '--viral-gene-required',
    default=False,
    is_flag=True,
    show_default=True,
    help=('require viral genes annotated, removing putative viral seqs '
            'with no genes annotated; this can reduce false positives at '
            'reasonable cost of sensitivity'),
)
@click.option(
    '--viral-gene-enrich-off',
    default=False,
    is_flag=True,
    show_default=True,
    help=('turn off the requirement of more viral than cellular genes '
            'for calling full sequence viral; this is useful when only '
            'using VirSorter2 to produce DRAMv input with viral sequence '
            'identified from other tools, or those trimmed by checkV'),
)
@click.option(
    '--seqname-suffix-off',
    default=False,
    is_flag=True,
    show_default=True,
    help=('turn off adding suffix (||full, ||{i}_partial, ||lt2gene) '
          'to sequence names; note that this might cause partial seqs '
          'from the same contig to have the same name; this option is '
          'could be used when you are sure there is one partial sequence '
          'at max from each contig'),
)
@click.option(
    '--provirus-off',
    default=False,
    is_flag=True,
    show_default=True,
    help=('turn off extracting provirus after classifying full contigs; '
            'Togetehr with lower --max-orf-per-seq, can speed up a run '
            'significantly'),
)
@click.option(
    '--max-orf-per-seq',
    default=-1,
    type=int,
    show_default=True,
    help=('max # of orf used for computing taxonomic feature; this option '
            'can only be used in --provirus-off mode; if # of orf in a seq '
            'exceeds the max limit, it is sub-sampled to this # '
            'to reduce computation')
)
@click.option(
    '--tmpdir',
    default='iter-0',
    help='directory name for intermediate files',
)
@click.option(
    '--rm-tmpdir',
    is_flag=True,
    default=False,
    show_default=True,
    help=('remove intermediate file directory (--tmpdir); More than 100 '
            'intermediate files are created for each run, so this is '
            'recommended when 100s of samples are processed to avoid '
            'exceeding file # quota for user')
)
@click.option(
    '--verbose',
    is_flag=True,
    default=False,
    show_default=True,
    help='shows detailed rules output',
)
@click.option(
    '--profile',
    default=None,
    help='snakemake profile e.g. for cluster execution.',
)
@click.option(
    '-n',
    '--dryrun',
    is_flag=True,
    default=False,
    show_default=True,
    help='Check rules to run and files to produce',
)
@click.option(
    '--use-conda-off',
    is_flag=True,
    default=False,
    show_default=True,
    help=('Stop using the conda envs (vs2.yaml) that come with this '
            'package and use what are installed in current system; '
            'Only useful when you want to install dependencies on '
            'your own with your own prefer version; this option '
            'works with the development version'),
)
@click.argument(
    'snakemake_args', 
    nargs=-1, 
    type=click.UNPROCESSED, 
)
def run_workflow(workflow, working_dir, db_dir, seqfile,
                 include_groups, jobs,  min_score, hallmark_required, 
                 hallmark_required_on_short, viral_gene_required, 
                 provirus_off, max_orf_per_seq, min_length, 
                 prep_for_dramv, tmpdir, rm_tmpdir, verbose, profile, 
                 dryrun, use_conda_off, snakemake_args, label, 
                 keep_original_seq, high_confidence_only, exclude_lt2gene, 
                 seqname_suffix_off, viral_gene_enrich_off):
    ''' Runs the virsorter main function to classify viral sequences

    This includes 3 steps: 1) preprocess, 2) feature extraction, and 3)
    classify. By default ("all") all steps are executed. The "classify"
    only run the 3) classify step without previous steps that are
    computationally heavy, good for rerunning with different filtering
    options (--min-score, --high-confidence-only, --hallmark-required,
    --hallmark-required-on-short, --viral-gene-required, --exclude-lt2gene).
    Most snakemake arguments can be appended to the command for more 
    info see 'snakemake --help'.
    '''

    # hard coded, need to change all "iter-0" to Tmpdir in smk
    tmpdir = 'iter-0'

    os.makedirs(working_dir, exist_ok=True)
    config_f = os.path.join(working_dir,'config.yaml')
    layout_f = os.path.join(working_dir, '.virsorter-layout.json')
    expected_layout = {'controller': 'snakemake9', 'checkpoints': 'split-marker-v1'}
    if os.path.exists(layout_f):
        with open(layout_f) as handle:
            if json.load(handle) != expected_layout:
                raise click.ClickException('Incompatible fork checkpoint layout; use a new work directory')
    elif os.path.exists(config_f):
        raise click.ClickException('Existing run has no compatible fork layout receipt; use a new work directory, not an in-place upgrade')
    else:
        temporary = layout_f + '.tmp'
        with open(temporary, 'w') as handle:
            json.dump(expected_layout, handle)
            handle.write('\n')
        os.replace(temporary, layout_f)

    if min_score > 1 or min_score < 0:
        logging.critical('--min-score needs to be between 0 and 1')
        sys.exit(1)
    if min_length < 0:
        logging.critical('--min-length needs to be >= 0')
        sys.exit(1)
    if jobs < 1:
        logging.critical('--jobs needs to be >= 1')
        sys.exit(1)

    if provirus_off:
        provirus = False
        if max_orf_per_seq != -1 and prep_for_dramv:
            mes = ('--max-orf-per-seq CAN NOT be used with '
                    '--prep-for-dramv; '
                    'outputs with ORFs subsampled are NOT '
                    'compatible with DRAMv')
    else:
        provirus = True
        max_orf_per_seq = -1

    if workflow == 'classify':
        if not os.path.exists(config_f):
            mes = 'No config.yaml dectected from previous run'
            logging.critical(mes)
            sys.exit(1)

        config = load_configfile(config_f)
        min_length_prev = config['MIN_LENGTH']
        if min_length != min_length_prev:
            mes = (
                '--min-length has changed from '
                f'{min_length_prev} to {min_length}; '
                'but --min-length has not effect on classify step; '
                'The whole pipeline has to be rerun if --min-length changes'
            )
            logging.critical(mes)
            sys.exit(1)

        if provirus != config['PROVIRUS']:
            mes = (
                '--provirus-off setting change found; '
                'The whole pipeline has to be rerun if --provirus-off changes'
            )
            logging.critical(mes)
            sys.exit(1)

        target_f = '{working_dir}/{tmpdir}/reclassify.trigger'.format(
                working_dir=working_dir,
                tmpdir=tmpdir,
        )
        try:
            subprocess.run(['touch', target_f], check=True)
        except subprocess.CalledProcessError as e:
            # removes the traceback
            #logging.critical(e)
            sys.exit(1)

        os.rename(config_f, os.path.join(working_dir,'config.yaml.bak'))

    make_config(
        db_dir=db_dir, seqfile=seqfile, include_groups=include_groups,
        threads=jobs, config_f=config_f, provirus=provirus,
        hallmark_required=hallmark_required,
        hallmark_required_on_short=hallmark_required_on_short,
        viral_gene_required=viral_gene_required,
        prep_for_dramv=prep_for_dramv,
        max_orf_per_seq=max_orf_per_seq, 
        tmpdir=tmpdir, min_length=min_length, min_score=min_score, 
        label=label, keep_original_seq=keep_original_seq,
        high_confidence_only=high_confidence_only, 
        exclude_lt2gene=exclude_lt2gene,
        seqname_suffix_off=seqname_suffix_off,
        viral_gene_enrich_off=viral_gene_enrich_off,
    )
    config = load_configfile(config_f)

    if db_dir == None:
        db_dir = config['DBDIR']

    cmd = snakemake_command(
        get_snakefile(), working_dir, jobs, configfile=config_f,
        conda_prefix=os.path.join(db_dir, 'conda_envs'),
        use_conda=not use_conda_off, profile=profile, dryrun=dryrun,
        verbose=verbose, targets=[workflow],
        force=[workflow] if workflow != 'all' else [], extra=snakemake_args,
    )
    try:
        execute_snakemake(cmd)
    except subprocess.CalledProcessError as e:
        logging.critical('Snakemake failed with exit status %s', e.returncode)
        sys.exit(e.returncode)

    # Retain native provenance/locks, including after dry-runs.
    if rm_tmpdir and not dryrun:
        to_remove = [tmpdir]
        for di in to_remove:
            _path = os.path.join(working_dir, di)
            shutil.rmtree(_path, ignore_errors=True)

# initialize
@cli.command(
    'setup',
    context_settings=dict(ignore_unknown_options=True),
    short_help='download reference files (~10GB) and install dependencies',
)
@click.option('-d',
    '--db-dir',
    help='diretory path for databases',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='number of simultaneous downloads',
)
@click.option(
    '-s',
    '--skip-deps-install',
    is_flag=True,
    show_default=True,
    help=('skip dependency installation (if you want to install on your '
            'own as in development version)'),
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def run_setup(db_dir,jobs, skip_deps_install, snakemake_args):
    '''Setup databases and install dependencies.
    
    Executes a snakemake workflow to download reference database files
    and validate based on their MD5 checksum, and install dependencies
    '''
    db_dir = os.path.abspath(db_dir)
    def setup_command(snakefile):
        return snakemake_command(
            get_snakefile(snakefile), db_dir, jobs,
            config={'Skip_deps_install': skip_deps_install},
            conda_prefix=os.path.join(db_dir, 'conda_envs'), extra=snakemake_args,
        )
    logging.info('Setting up VirSorter2 database; this only needs to be done once.')
    try:
        execute_snakemake(setup_command('rules/setup.smk'))
    except subprocess.CalledProcessError:
        logging.info('First download attempt failed; trying the fallback source.')
        try:
            execute_snakemake(setup_command('rules/setup-retry.smk'))
        except subprocess.CalledProcessError as e:
            logging.critical(e)
            sys.exit(e.returncode)

# train feature
@cli.command(
    'train-feature',
    context_settings=dict(ignore_unknown_options=True),
    short_help='subcommand for training feature of customized classifier',
)
@click.option('-w',
    '--working-dir',
    help='output directory',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    '--seqfile',
    help=('genome sequence file for training; for file pattern globbing, '
            'put quotes around the pattern, eg. "fasta-dir/*.fa"'),
    type=str,
    required=True,
    multiple=True,
)
@click.option(
    '--hmm',
    help=('customized viral HMMs for training; default to the one used in '
            'VirSorter2'),
    type=click.Path(resolve_path=True),
)
@click.option(
    '--hallmark',
    type=click.Path(resolve_path=True),
    help=('hallmark gene hmm list from -hmm for training (a tab separated '
            'file with three columns: 1. hmm name, 2. gene name of hmm, '
            '3. hmm bit score cutoff); default to the one used for '
            'dsDNAphage in VirSorter2'),
)
@click.option(
    '--prodigal-train',
    help=('customized training db from prodigal; default to the one used '
            'in prodigal --meta mode'),
    type=click.Path(resolve_path=True),
)
@click.option(
    '--frags-per-genome',
    default=5,
    type=int,
    show_default=True,
    help='number of random DNA fragments collected from each genome',
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='max # of jobs in parallel',
)
@click.option(
    '--min-length',
    default=1000,
    type=int,
    show_default=True,
    help='minimum size of random DNA fragment for training',
)
@click.option(
    '--max-orf-per-seq',
    default=20,
    type=int,
    show_default=True,
    help=('Max # of orf used for computing taxonomic features; if # of orf '
            'in a seq exceeds the max limit, it is sub-sampled to this # '
            'to reduce computation; to turn off this, set it to -1'),
)
@click.option(
    '--genome-as-bin',
    default=False,
    is_flag=True,
    show_default=True,
    help=('if applied, each file (genome bin) is a genome in --seqfile, '
            'else each sequence is a genome'),
)
@click.option(
    '--use-conda-off',
    is_flag=True,
    default=False,
    show_default=True,
    help=('Stop using the conda envs (vs2.yaml) that come with this package '
            'and use what are installed in current system; Only useful when '
            'you want to install dependencies on your own with your own '
            'preferred versions'),
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def train_feature(working_dir, seqfile, hmm, hallmark, prodigal_train, 
                    frags_per_genome, min_length, max_orf_per_seq, 
                    genome_as_bin, jobs, use_conda_off, snakemake_args):
    '''Training features for customized classifier.
    
    Executes a snakemake workflow to do the following:
    1) prepare random DNA fragments from viral and nonviral genome data 
    2) extract feature from random DNA fragments to make ftrfile
    '''

    DEFAULT_CONFIG = get_default_config()

    cwd = os.getcwd()
    lis = []
    pat_lis = []
    for pat in seqfile:
        # only works in linux
        if pat.startswith('/'):
            new_pat = pat
        else:
            new_pat = '{}/{}'.format(cwd, pat)
        fs = glob.glob(pat)
        lis.extend(fs)
        pat_lis.append(new_pat)
    
    if len(lis) == 0:
        mes = 'No files match {}'.format(viral_seqfile)
        logging.critical(mes)
        sys.exit(1)
    else:
        mes = '{} seqfiles are used for training features'.format(len(lis))
        logging.info(mes)

    if hmm == None:
        hmm = 'NA'
    if hallmark == None:
        hallmark = 'NA'

    if prodigal_train == None:
        prodigal_train = 'NA'

    cmd = snakemake_command(
        get_snakefile('rules/train-feature.smk'), working_dir, jobs,
        config={'Viral_seqfile': ' '.join(pat_lis), 'Hmm': hmm,
                'Hallmark': hallmark, 'Rbs': prodigal_train,
                'Min_length': min_length, 'Max_orf_per_seq': max_orf_per_seq,
                'Viral_genome_as_bin': genome_as_bin, 'Fragments_per_genome': frags_per_genome},
        conda_prefix=os.path.join(DEFAULT_CONFIG['DBDIR'], 'conda_envs'),
        use_conda=not use_conda_off, extra=snakemake_args,
    )
    try:
        execute_snakemake(cmd)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)

# train model
@cli.command(
    'train-model',
    context_settings=dict(ignore_unknown_options=True),
    short_help='subcommand for training customized classifier model',
)
@click.option('-w',
    '--working-dir',
    help='output directory',
    type=click.Path(dir_okay=True,writable=True,resolve_path=True),
    required=True
)
@click.option(
    '--viral-ftrfile',
    help='viral genome feature file for training',
    type=click.Path(resolve_path=True),
    required=True,
)
@click.option(
    '--nonviral-ftrfile',
    help='nonviral genome feature file for training',
    type=click.Path(resolve_path=True),
    required=True,
)
@click.option(
    '-j',
    '--jobs',
    default=multiprocessing.cpu_count(),
    type=int,
    show_default=True,
    help='number of threads for classier',
)
@click.option(
    '--balanced',
    is_flag=True,
    default=False,
    show_default=True,
    type=bool,
    help='random undersample the larger to the size of the smaller feature file'
)
@click.option(
    '--use-conda-off',
    is_flag=True,
    default=False,
    show_default=True,
    help=('Stop using the conda envs (vs2.yaml) that come with this package '
            'and use what are installed in current system; Only useful when '
            'you want to install dependencies on your own with your own '
            'prefer versions'),
)
@click.argument('snakemake_args', nargs=-1, type=click.UNPROCESSED)
def train_model(working_dir, viral_ftrfile, nonviral_ftrfile, balanced, 
                jobs, use_conda_off, snakemake_args):
    '''Training customized classifier model.
    '''

    DEFAULT_CONFIG = get_default_config()

    if balanced == None:
        balanced = False
    cmd = snakemake_command(
        get_snakefile('rules/train-model.smk'), working_dir, jobs,
        config={'Viral_ftrfile': viral_ftrfile, 'Nonviral_ftrfile': nonviral_ftrfile,
                'Balanced': balanced, 'Jobs': jobs},
        conda_prefix=os.path.join(DEFAULT_CONFIG['DBDIR'], 'conda_envs'),
        use_conda=not use_conda_off, extra=snakemake_args,
    )
    try:
        execute_snakemake(cmd)
    except subprocess.CalledProcessError as e:
        sys.exit(e.returncode)

# config management
@cli.command(
    'config',
    context_settings=dict(ignore_unknown_options=True),
    #no_args_is_help=True,
    short_help='subcommand for configuration management',
)
@click.option(
    '--show',
    is_flag=True,
    default=False,
    show_default=True,
    help='show all configuration values',
)
@click.option(
    '--show-source',
    is_flag=True,
    default=False,
    show_default=True,
    help='show path of the configuration file',
)
@click.option(
    '--init-source',
    is_flag=True,
    default=False,
    show_default=True,
    help='initialize configuration file',
)
@click.option(
    '--db-dir',
    type=click.Path(),
    help='directory for databases; required for --init-source',
)
@click.option(
    '--set',
    help=('set KEY to VAL with the format: KEY=VAL; for nested dict in YAML '
            'use KEY1.KEY2=VAL (e.g. virsorter config --set '
            'GROUP_INFO.RNA.MIN_GENOME_SIZE=2000)'),
)
@click.option(
    '--get',
    help=('the value of a KEY (e.g. virsorter config '
            '--get GROUP_INFO.RNA.MIN_GENOME_SIZE'),
)
def config(show, show_source, init_source, db_dir, set, get):
    '''CLI for managing configurations.

    There are many configurations kept in "template-config.yaml" in source 
    code directory or "~/.virsorter" (when source code directory is not 
    writable for user). This file can located with 
    `virsorter config --show-source`. You can set the configurations with 
    `virsorter config --set KEY=VAL`. Alternative, you can edit in the 
    configuration file ("template-config.yaml") directly.
    '''

    from virsorter.config import (
            TEMPLATE, SRC_CONFIG_DIR, 
            USER_CONFIG_DIR, init_config_template
    )

    if init_source:
        if db_dir == None:
            mes = '--db-dir is required for --init-source'
            logging.critical(mes)
            sys.exit(1)
        else:
            if not os.path.isdir(db_dir):
                mes = (f'--db-dir {db_dir} does NOT exist yet; Make sure it '
                        'is created later\n')
                logging.warning(mes)

            db_dir = os.path.abspath(db_dir)
            init_config_template(SRC_CONFIG_DIR, USER_CONFIG_DIR, db_dir)
            sys.exit(0)

    if not os.path.isfile(TEMPLATE):
        mes = ('config file "template-config.yaml" has not been '
                'initialized yet; Please use '
                '`virsorter config --init-source --db-dir PATH` to initialize')
        logging.critical(mes)
        sys.exit(1)

    config = get_default_config()

    if show:
        YAML().dump(config, sys.stdout)
        sys.exit(0)

    if show_source:
        mes = f'config file path: {TEMPLATE}\n'
        sys.stdout.write(mes)
        sys.exit(0)

    if get != None:
        s = get
        lis = [var.strip() for var in s.split(',')]
        for var in lis:
            temp = config
            for i in var.split('.'):
                i = i.strip()
                try:
                    temp = temp[i]
                except KeyError as e:
                    mes = f'{i} is not a key in config file ({TEMPLATE})'
                    logging.critical(mes)
                    sys.exit(1)

            mes = f'{var}: {temp}\n'
            sys.stdout.write(mes)

        sys.exit(0)

    if set != None:
        s = set
        lis = [item.strip() for item in s.split(',')]
        for item in lis:
            temp = config
            var, val = item.split('=')
            var = var.strip()
            val = val.strip()
            keys = [key.strip() for key in var.split('.')]
            for i in range(len(keys)):
                if i == (len(keys) - 1):
                    # stop at 2nd last key
                    break
                key = keys[i]
                try:
                    temp = temp[key]
                except KeyError as e:
                    mes = f'{key} is not a key in config file ({TEMPLATE})'
                    logging.critical(mes)
                    sys.exit(1)

            last_key = keys[-1]

            try:
                old_val = temp[last_key]
                if isinstance(old_val, int):
                    try:
                        val = int(val)
                    except ValueError as e:
                        mes = f'{var} is supposed to be an integer'
                        logging.critical(mes)
                        sys.exit(1)
                elif isinstance(old_val, float):
                    val = float(val)
                    try:
                        val = float(val)
                    except ValueError as e:
                        mes = f'{var} is supposed to be a float'
                        logging.critical(mes)
                        sys.exit(1)
                # only convert to abspath when the old one exists
                #   since sometimes just want to set relative path
                elif os.path.exists(old_val):
                    val = os.path.abspath(val)
                temp[last_key] = val
            except KeyError as e:
                mes = f'{last_key} is not a key in config file ({TEMPLATE})'
                logging.critical(mes)
                sys.exit(1)

            mes = f'{var}: {old_val} ==> {val}\n'
            sys.stdout.write(mes)
            with open(TEMPLATE, 'w') as fw:
                YAML().dump(config, fw)

        sys.exit(0)
