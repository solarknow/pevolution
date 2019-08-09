import os

# command constants
ALIGN = 'clustalw -align -infile={input}.edit -outfile={output}.ed -output={fmt} -quiet'
PROTTEST = 'java -jar prottest/prottest-3.4.2.jar -i {infile} -o {outfile}' + \
               '-all-distributions -all -S 1 -threads {num_procs} -BIC'


def clustal_align(infile, outfile, fmt='nexus'):
    formatted = ALIGN.format(input=infile, output=outfile, fmt=fmt)
    os.system(formatted)


def run_prottest(infile, outfile, procs):
    os.system(PROTTEST.format(infile=infile, outfile=outfile, num_procs=procs))