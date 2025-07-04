import os
from dataclasses import dataclass


@dataclass
class PathDef:
    filename: str
    dir: str

    def __str__(self):
        os.makedirs(self.dir, exist_ok=True)
        return self.dir + os.sep + self.filename


@dataclass
class OrthosPath(PathDef):
    dir: str = 'Orthos'


@dataclass
class XMLPath(PathDef):
    dir: str = 'XML'


@dataclass
class DictsPath(PathDef):
    dir: str = 'dicts'


@dataclass
class DataPath(PathDef):
    dir: str = 'Data'


@dataclass
class AlignsPath(PathDef):
    dir: str = 'aligns'


@dataclass
class ProtPath(PathDef):
    dir: str = 'Prot'


@dataclass
class BayesPath(PathDef):
    dir: str = 'Bayes'


@dataclass
class MLPath(PathDef):
    dir: str = 'ML'


@dataclass
class ReportsPath(PathDef):
    dir: str = 'Reports'


def resolve_prottest_path():
    if 'prottest' in os.listdir():
        listdir = filter(lambda x: x[-3:] == 'jar', os.listdir('prottest'))
        return PathDef(list(listdir)[0], 'prottest')
    return None


@dataclass
class BlastThresholds:
    arch: float
    bac: float
    euk: float

class EVals:
    CLOSE = 1e-10
    SIMILAR = 1e-5
    DISTANT = 5

class Domains:
    EUKARYOTA="Eukaryota"
    BACTERIA = "Bacteria"
    ARCHAEA = "Archaea"

def organism_list(domain):
    arch_list = ['Haloferax volcanii', 'Sulfolobus tokodaii', 'Methanococcus aeolicus',
                 'Methanobrevibacter smithii', 'Thermococcus sibiricus', 'Archaeoglobus fulgidus',
                 'Nanoarchaeum equitans', 'Thermoplasma acidophilum']
    bac_list = ['Bacillus subtilis', 'Escherichia coli', 'Gemmata obscuriglobus',
                'Rickettsia prowazekii', 'Agrobacterium tumefaciens', 'Prosthecobacter dejongeii',
                'Anabaena variabilis', 'Thermotoga maritima', 'Verrucomicrobium spinosum']
    euk_list = ['Homo sapiens', 'Drosophila melanogaster', 'Oryza sativa',
                'Trypanosoma brucei', 'Plasmodium falciparum', 'Saccharomyces cerevisiae',
                'Neurospora crassa', 'Arabidopsis thaliana']
    if domain == Domains.EUKARYOTA:
        return euk_list
    elif domain == Domains.BACTERIA:
        return bac_list
    elif domain == Domains.ARCHAEA:
        return arch_list
    else:
        return arch_list + bac_list + euk_list

def domain_thresholds(domain):
    if domain == Domains.EUKARYOTA:
        return BlastThresholds(arch=EVals.DISTANT, bac=EVals.DISTANT, euk=EVals.CLOSE)
    elif domain == Domains.ARCHAEA:
        return BlastThresholds(arch=EVals.CLOSE, bac=EVals.SIMILAR, euk=EVals.DISTANT)
    elif domain == Domains.BACTERIA:
        return BlastThresholds(arch=EVals.SIMILAR, bac=EVals.CLOSE, euk=EVals.DISTANT)
    else:
        return BlastThresholds(arch=EVals.DISTANT, bac=EVals.DISTANT, euk=EVals.DISTANT)


NEXUS = '\n'.join(['#NEXUS', 'begin data;',
                   'dimensions ntax={num_taxa} nchar={num_char};',
                   'format datatype={type} interleave=no gap=-;',
                   'matrix', '', ''])


def nexus_fmt(num_seq, seq_len, data_type='protein'):
    return NEXUS.format(num_taxa=num_seq, num_char=seq_len, type=data_type)

