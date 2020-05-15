import os

from Bio import Entrez


def set_email(email):
    Entrez.email = email


def fetch_protein(acc, fmt='fasta'):
    """
    returns a file-like handle of a Entrez search for acc accession number and returns in fmt format
    """
    fetched = Entrez.efetch(db='protein', id=str(acc), rettype=fmt)
    return fetched.readlines()


def fetch_fasta(acc):
    """
    prints the sequence in fasta format to file
    """
    hand = fetch_protein(acc)
    org = fetch_organism(acc)
    string = '>{spec}: {domain} {id}\n'.format(
        spec=org[0],
        domain=org[1],
        id=acc
    )
    for line in hand[1:]:
        string += line
        if line == '\n':
            os.makedirs('Orthos', exist_ok=True)
            new_file = 'Orthos' + os.sep + acc + '.fasta'
            with open(new_file, 'w') as writ_file:
                writ_file.write(string)
            return new_file


def fetch_definition(acc):
    """
    Fetchs definition of specified gene product
    """
    hand = fetch_protein(acc, 'gp')
    record = '  '.join(line.strip() for line in hand)
    rec = record.split('  ')
    recun = [s for s in rec if s != '']
    i = 0
    name = ''
    while i < len(recun):
        if recun[i].endswith('DEFINITION'):
            while i:
                i += 1
                if recun[i].startswith('ACCESSION'):
                    return name.strip()
                name += recun[i] + ' '
        i += 1


def fetch_organism(acc):
    """
    Fetchs the organism from which acc accession number of a protein came from
    """
    retu = []
    hand = fetch_protein(acc, 'gp')
    record = '  '.join(line.strip() for line in hand)
    rec = record.split('  ')
    recun = [s for s in rec if s != '']
    i = 0
    while i < len(recun):
        if recun[i] == 'DEFINITION':
            multisp = recun[i + 1].startswith("MULTISPECIES:")
            if not multisp:
                tab = recun[i + 1].strip().split()
                if not recun[i + 2] == 'ACCESSION':
                    tab += recun[i + 2].strip().split()
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
                return [recun[i + 11] + ' multispecies', recun[i + 14].split(';')[0]]

        if recun[i] == 'ORGANISM':
            try:
                ret = [recun[i + 1].split()[0] + ' ' + recun[i + 1].split()[1]]
            except IndexError:
                return retu
            if ret[0].endswith('.'):
                ret[0] += ' ' + recun[i + 1].split()[2]
            ret.append(recun[i + 2].split(';')[0])
            return ret

        i += 1
