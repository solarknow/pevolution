import os
from dataclasses import dataclass


@dataclass
class PathDef:
    filename: str
    dir: str
    os.makedirs(dir, exist_ok=True)

    def __str__(self):
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
