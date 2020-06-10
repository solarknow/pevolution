import os
import subprocess

from Bio import Entrez

from constants import ORTHOS_PATH


class Protein:
    def __init__(self, gp):
        """
        Taking the lines of a GenPept file, extracts and sets class fields
        :param gp: list of strings, lines from a GenPept file
        """
        self.seq = ''
        self.definition = ''
        self.domain = ''
        self.binomial = ''
        self.accession = ''
        processed_gp = [s for s in '  '.join(line.strip() for line in gp).split('  ') if s != '']
        # print(processed_gp)
        for i, word in enumerate(processed_gp):
            if word == 'DEFINITION':
                for sub_word in processed_gp[i + 1:]:
                    if sub_word == 'ACCESSION':
                        break
                    else:
                        self.definition += sub_word
            elif word == 'ACCESSION':
                self.accession = processed_gp[i + 1].strip()
            elif word == 'ORGANISM':
                self.binomial = processed_gp[i + 1]
                self.domain = processed_gp[i + 2].split(';')[0]
            elif word == 'ORIGIN':
                for sub_word in processed_gp[i + 1:-1]:
                    split_word = sub_word.split()
                    self.seq += ''.join(split_word[1:]).upper() + '\n'

    @classmethod
    def manual(cls, seq, definition, domain, binomial, accession):
        obj = cls.__new__(cls)
        super(Protein, obj).__init__()
        obj.seq = seq
        obj.definition = definition
        obj.domain = domain
        obj.binomial = binomial
        obj.accession = accession
        return obj

    def __str__(self):
        return '''Accession: {acc}
        Definition: {definition}
        Domain: {dom}
        Binomial: {bin}
        Sequence: {seq}'''.format(acc=self.accession, definition=self.definition, dom=self.domain, bin=self.binomial,
                                  seq=self.seq)

    def __eq__(self, other):
        return self.seq == other.seq and \
               self.accession == other.accession and \
               self.definition == other.definition and \
               self.domain == other.domain and \
               self.binomial == other.binomial


def set_email(email):
    Entrez.email = email


def fetch_protein(accession):
    """
    Creates a Protein object of a Entrez search for acc accession number
    :param accession: Accession ID or GI number
    :return: Protein instance of fetched data for Accession
    """
    fetched = Entrez.efetch(db='protein', id=str(accession), rettype='gp')
    return Protein(fetched.readlines())


def write_fasta(protein):
    """
    prints the sequence in fasta format to file
    :param protein: Protein instance to write to file
    :return: path to file that was written
    """
    string = '>{spec}: {domain} {id}\n{seq}'.format(
        spec=protein.binomial,
        domain=protein.domain,
        id=protein.accession,
        seq=protein.seq
    )
    os.makedirs(ORTHOS_PATH, exist_ok=True)
    new_file = ORTHOS_PATH + protein.accession + '.fasta'
    if not os.path.exists(new_file):
        with open(new_file, 'w') as write_file:
            write_file.write(string)
    return new_file


def remote_blast(query, threshold, outfile, organism):
    """
    Runs a remote BLAST query.
    :param query: File path to the sequence to use as a query
    :param threshold: E-value threshold
    :param outfile: File path to output file
    :param organism: Binomial name of organism to query against
    :return: None
    """
    blast_query = ['blastp', '-db', 'nr', '-query', query, '-evalue', str(threshold), '-out',
                   outfile, '-outfmt', str(5), '-entrez_query', organism + '[ORGN]', '-remote']
    return subprocess.run(blast_query)
