#!/usr/bin/env python

try:
  from setuptools import setup
except ImportError:
  from distutils.core import setup
from urllib.request import urlretrieve
import tarfile, zipfile


setup(name='pevolution',
      version='1.0',
      description='This is a pipeline to find, align, and find trees for putatively related proteins.',
      author='Mihir Sarwade',
      author_email='mihir.sarwade@effem.com',
#      packages=['pevolution'],
      install_requires=['psutil >= 1.2.1','biopython >= 1.64'])

urlretrieve("https://prank-msa.googlecode.com/files/prank.linux64.140110.tgz",'/tmp/prank.tgz')
urlretrieve("http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip",'/tmp/PhyML.zip')
tarfile.open('/tmp/prank.tgz').extractall(path='~/prank')
zipfile.open('/tmp/PhyML').extractall(path='~/phyml')
tarfile.open('prottest*.tar.gz').extractall(path='~/prottest')