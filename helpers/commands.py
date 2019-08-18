import os

# command constants
ALIGN = 'clustalw -align -infile={input}.edit -outfile={output}.ed -output={fmt} -quiet'
PROTTEST = 'java -jar prottest' + os.sep + 'prottest-3.4.2.jar -i {infile} -o {outfile}' + \
           '-all-distributions -all -S 1 -threads {num_procs} -BIC'
BLAST = 'blastp -db nr -query Orthos' + os.sep + '{seed}.fasta -evalue {threshold} -out XML' + \
        os.sep + '{xml_file}.xml -outfmt 5 -entrez_query \"{org}[ORGN]\" -use_sw_tback -remote'


def clustal_align(infile, outfile, fmt='nexus'):
    formatted = ALIGN.format(input=infile, output=outfile, fmt=fmt)
    os.system(formatted)


def run_prottest(infile, outfile, procs):
    os.system(PROTTEST.format(infile=infile, outfile=outfile, num_procs=procs))


def run_blast(seed, thresh, dum, org):
    os.system(BLAST.format(seed=seed, threshold=thresh, xml_file=dum, org=org))
