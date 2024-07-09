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
class OrthoPath(PathDef):
    dir: str = 'Orthos'


@dataclass
class XMLPath(PathDef):
    dir: str = 'XML'


@dataclass
class DictsPath(PathDef):
    dir: str = 'dicts'


def resolve_prottest_path():
    if 'prottest' in os.listdir():
        listdir = filter(lambda x: x[-3:] == 'jar', os.listdir('prottest'))
        return PathDef(listdir[0], 'prottest')
