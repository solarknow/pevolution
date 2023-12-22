import os

from Bio import Entrez


def set_email(email):
    Entrez.email = email


def fetch_protein(acc, fmt='fasta'):
    """returns a file-like handle of an Entrez search for acc accession number and returns in fmt format"""
    fetched = Entrez.efetch(db='protein', id=str(acc), rettype=fmt)
    return fetched.readlines()


def fetch_fasta(acc):
    """prints the sequence in fasta format to file"""
    hand = fetch_protein(acc)
    string = ''
    org = fetch_organism(acc)
    string += '>' + org[0] + ': ' + org[1] + ' ' + acc + '\n'
    for line in hand[1:]:
        string += line
        if line == '\n':
            os.makedirs('Orthos', exist_ok=True)
            with open('Orthos' + os.sep + acc + '.fasta', 'w') as writ_file:
                writ_file.write(string)
            return 'Orthos' + os.sep + acc + '.fasta'


def fetch_definition(acc):
    """Fetches definition of specified gene product"""
    protein = fetch_protein(acc, 'gp')
    record = ''
    for line in protein:
        record += line.strip() + '  '
    rec = record.split('  ')
    non_null_recs = [s for s in rec if s]
    i = 0
    name = ''
    while i < len(non_null_recs):
        if non_null_recs[i].endswith('DEFINITION'):
            while i:
                i += 1
                if non_null_recs[i].startswith('ACCESSION'):
                    return name.strip()
                name += non_null_recs[i] + ' '
        i += 1


def fetch_organism(acc):
    """Fetches the organism from which acc accession number of a protein came from"""
    protein = fetch_protein(acc, 'gp')
    record = ''
    for line in protein:
        record += line.strip() + '  '
    rec = record.split('  ')
    non_null_rec = [s for s in rec if s]
    i = 0
    ret_list = []
    retu = ['']
    while i < len(non_null_rec):
        if non_null_rec[i] == 'DEFINITION':
            is_multispecies = non_null_rec[i + 1].startswith("MULTISPECIES:")
            if not is_multispecies:
                tab = non_null_rec[i + 1].strip().split()
                if not non_null_rec[i + 2] == 'ACCESSION':
                    tab += non_null_rec[i + 2].strip().split()
                for k in range(len(tab)):
                    if tab[k].startswith('[') and len(tab[k]) > 5:
                        retu = [tab[k][1:]]
                        while k < len(tab):
                            k += 1
                            if tab[k].endswith('].'):
                                retu[0] += ' ' + tab[k][0:tab[k].index(']')]
                                retu.append('Unknown.')
                                break
                            else:
                                retu[0] += ' ' + tab[k]
            else:
                return [non_null_rec[i + 11] + ' multispecies', non_null_rec[i + 14].split(';')[0]]

        if non_null_rec[i] == 'ORGANISM':
            try:
                ret_list.append(non_null_rec[i + 1].split()[0] + ' ' + non_null_rec[i + 1].split()[1])
            except IndexError:
                return retu
            if ret_list[0].endswith('.'):
                ret_list[0] += ' ' + non_null_rec[i + 1].split()[2]
            ret_list.append(non_null_rec[i + 2].split(';')[0])
            return ret_list

        i += 1
