import subprocess
from helpers.constants import OrthosPath, XMLPath, resolve_prottest_path

# command constants
ALIGN = ['clustalw', '-align', '-infile={input}.edit', '-outfile={output}.ed', '-output={fmt}', '-quiet']
PROTTEST = ['java', '-jar', str(resolve_prottest_path()), '-i', '{infile}', '-o', '{outfile}',
            '-all-distributions', '-all', '-S', 1, '-threads', '{num_procs}', '-BIC']
BLAST = ['blastp', '-db', 'nr', '-query', str(OrthosPath('{seed}.fasta')), '-evalue', '{threshold}',
         '-out', str(XMLPath('{xml_file}.xml')), '-outfmt', '5', '-entrez_query', '\"{org}[ORGN]\"',
         '-use_sw_tback', '-remote']
BAYES = ['mb', '{str(cmdfile)}']
PRANK = ['prank', '-d={infile}', '-o={outfile}', '-f=nexus', '-quiet']
PHYML = ['phyml', '-i', '{str(infile)}', '-d', 'aa', '-b', 100, '-m', '{model}',
         '-f', 'e', '-s', 'BEST', '-u', '{str(outfile)}', '-o', 'tl']
DOM = ['python3', 'dom.py', '{out}', '{query}', '{dom}', '{phyml}']


def format_run(cmd, **kwargs):
    subprocess.call([c.format(**kwargs) for c in cmd])


def clustal_align(infile, outfile, fmt='nexus'):
    format_run(ALIGN, input=infile, output=outfile, fmt=fmt)


def run_prottest(infile, outfile, procs):
    format_run(PROTTEST, infile=infile, outfile=outfile, num_procs=procs)


def run_blast(seed, thresh, dum, org):
    format_run(BLAST, seed=seed, threshold=thresh, xml_file=dum, org=org)


def run_prank(infile, outfile):
    format_run(PRANK, infile=infile, outfile=outfile)


def run_phyml(infile, model, outfile, v=None, a=None, ori=None):
    CMD = PHYML
    if v:
        CMD += ['-v', v]
    if a:
        CMD += ['-a', a]
    if ori:
        CMD.append('-ori')
    format_run(CMD, infile=infile, outfile=outfile, model=model)


def run_dom_file(out, query, dom, phyml=False):
    format_run(DOM, out=out, query=query, dom=dom, phyml=phyml)


def run_bayes(infile):
    format_run(BAYES, cmdfile=infile)
