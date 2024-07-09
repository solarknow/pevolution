import subprocess
from helpers.constants import OrthoPath, XMLPath, resolve_prottest_path

# command constants
ALIGN = ['clustalw', '-align', '-infile={input}.edit', '-outfile={output}.ed', '-output={fmt}', '-quiet']
PROTTEST = ['java', '-jar', str(resolve_prottest_path()), '-i', '{infile}', '-o', '{outfile}',
            '-all-distributions', '-all', '-S', 1, '-threads', '{num_procs}', '-BIC']
BLAST = ['blastp', '-db', 'nr', '-query', str(OrthoPath('{seed}.fasta')), '-evalue', '{threshold}',
         '-out', str(XMLPath('{xml_file}.xml')), '-outfmt', '5', '-entrez_query', '\"{org}[ORGN]\"',
         '-use_sw_tback', '-remote']


def format_run(cmd, **kwargs):
    subprocess.call([c.format(**kwargs) for c in cmd])


def clustal_align(infile, outfile, fmt='nexus'):
    format_run(ALIGN, input=infile, output=outfile, fmt=fmt)


def run_prottest(infile, outfile, procs):
    format_run(PROTTEST, infile=infile, outfile=outfile, num_procs=procs)


def run_blast(seed, thresh, dum, org):
    format_run(BLAST, seed=seed, threshold=thresh, xml_file=dum, org=org)
