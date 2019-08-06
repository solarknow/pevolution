import os

from Bio import Entrez


def set_email(email):
    Entrez.email = email


def fetch_protein(acc, fmt='fasta'):
    """returns a file-like handle of a Entrez search for acc accession number and returns in fmt format"""
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
            if not os.path.exists('Orthos'):
                os.mkdir('Orthos')
            with open('Orthos/' + acc + '.fasta', 'w') as writ_file:
                writ_file.write(string)
            return 'Orthos/' + acc + '.fasta'


def fetch_definition(acc):
    """Fetches definition of specified gene product"""
    hand = fetch_protein(acc, 'gp')
    records = '  '.join([line.strip() for line in hand])
    filtered_records = [record for record in records.split('  ') if record]
    i = 0
    name = ''
    while i < len(filtered_records):
        if filtered_records[i].endswith('DEFINITION'):
            while i:
                i += 1
                if filtered_records[i].startswith('ACCESSION'):
                    return name.strip()
                name += filtered_records[i] + ' '
        i += 1


def fetch_organism(acc):
    """Fetches the organism from which acc accession number of a protein came from"""
    hand = fetch_protein(acc, 'gp')
    records = '  '.join([line.strip() for line in hand])
    filtered_records = [record for record in records.split('  ') if record]
    i = 0
    retu = ['']
    while i < len(filtered_records):
        if filtered_records[i] == 'DEFINITION':
            multispecies = filtered_records[i + 1].startswith("MULTISPECIES:")
            if not multispecies:
                tab = filtered_records[i + 1].strip().split()
                if not filtered_records[i + 2] == 'ACCESSION':
                    tab += filtered_records[i + 2].strip().split()
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
                return [filtered_records[i + 11] + ' multispecies', filtered_records[i + 14].split(';')[0]]

        if filtered_records[i] == 'ORGANISM':
            try:
                ret = [filtered_records[i + 1].split()[0] + ' ' + filtered_records[i + 1].split()[1]]
            except IndexError:
                return retu
            if ret[0].endswith('.'):
                ret[0] += ' ' + filtered_records[i + 1].split()[2]
            ret.append(filtered_records[i + 2].split(';')[0])
            return ret

        i += 1
