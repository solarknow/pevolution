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
        return PathDef(listdir[0], 'prottest')
