import os
import subprocess

from helpers.constants import OrthosPath, XMLPath, resolve_prottest_path

# command constants
ALIGN = [
    "clustalw",
    "-align",
    "-infile={input}.edit",
    "-outfile={output}.ed",
    "-output={fmt}",
    "-quiet",
]
PROTTEST = [
    "java",
    "-jar",
    str(resolve_prottest_path()),
    "-i",
    "{infile}",
    "-o",
    "{outfile}",
    "-all-distributions",
    "-all",
    "-S",
    "1",
    "-BIC",
]
BLAST = [
    "blastp",
    "-db",
    "{db}",
    "-query",
    str(OrthosPath("{seed}.fasta")),
    "-evalue",
    "{threshold}",
    "-out",
    str(XMLPath("{xml_file}.xml")),
    "-outfmt",
    5,
]
CREATE_DB = ["makeblastdb", "-in", "{infile}", "-dbtype", "prot", "-out", "{db}", "-parse_seqids"]
BAYES = ["mb", "{cmdfile}"]
PRANK = ["prank-msa" + os.sep + "prank", "-d={infile}", "-o={outfile}", "-f=nexus", "-showall", "-quiet"]
PHYML = [
    os.sep.join(["phyml", "src", "phyml"]),
    "-i",
    "{infile}",
    "-d",
    "aa",
    "-b",
    "100",
    "-m",
    "{model}",
    "-f",
    "e",
    "-u",
    "{outfile}",
    "-o",
    "tl",
]
DOM = ["python3", "dom.py", "{out}", "{query}", "{dom}", "{phyml}"]


def format_run(cmd, **kwargs):
    fmt = [str(c).format(**kwargs) for c in cmd]
    print(" ".join(fmt))
    subprocess.run(fmt, check=True)


def clustal_align(infile, outfile, fmt="nexus"):
    format_run(ALIGN, input=infile, output=outfile, fmt=fmt)


def run_prottest(infile, outfile):
    format_run(PROTTEST, infile=infile, outfile=outfile)


def run_blast(seed, thresh, dum, org, db="nr"):
    if db == "nr":
        new_blast = BLAST + ["-remote", "-entrez_query", f'"{org}[ORGN]"']
    else:
        new_blast = BLAST
    format_run(new_blast, seed=seed, threshold=thresh, xml_file=dum, org=org, db=db)


def run_prank(infile, outfile):
    format_run(PRANK, infile=infile, outfile=outfile)


def run_phyml(infile, model, outfile, v=None, a=None):
    cmd = list(PHYML)
    if v:
        cmd += ["-v", v]
    if a:
        cmd += ["-a", a]
    format_run(cmd, infile=infile, outfile=outfile, model=model)


def run_domain_file(out, query, dom, phyml=False):
    format_run(DOM, out=out, query=query, dom=dom, phyml=phyml)


def run_bayes(cmdfile):
    format_run(BAYES, cmdfile=cmdfile)
