from Bio.Blast import NCBIXML

from Utils import FetchUtil, SeqUtil
from helpers.constants import OrthosPath, DataPath


def xml_parse_and_extract_accession_numbers(xml_path: str) -> list[str]:
    """Returns a list of accession numbers that are contained within a Blast XML"""
    print("blasting done. XML file exists, parsing")
    with open(xml_path) as q1output:
        parse = NCBIXML.parse(q1output)
        acc = []
        for lin in parse:
            for align in lin.alignments:
                for hsp in sorted(align.hsps, key= lambda x: x.bit_score, reverse=True):
                    if (hsp.positives / float(hsp.align_length)) >= .4 and (
                            float(hsp.align_length) / len(hsp.query)) > .25:
                        acc.append(align.title.split('|')[1])
    return acc


def merge_domain_fastas(outfile, acc_dict):
    """Bring together FASTAs, defined by acc_dict"""
    for a in acc_dict.values():
        FetchUtil.fetch_fasta(a[0])
        original_fasta = OrthosPath(a[0] + '.fasta')
        modified_fasta = OrthosPath(a[0] + '.mod.fasta')
        with open(str(original_fasta)) as fil:
            fil_arr = fil.readlines()

        with open(str(modified_fasta), 'w') as fil:
            for i in range(len(fil_arr)):
                if i == 0:
                    fil.write(fil_arr[i].strip() + ' ' + a[1] + '  ' + a[2] + '\n')
                else:
                    fil.write(fil_arr[i])

        SeqUtil.append_sequences(DataPath(outfile), modified_fasta)