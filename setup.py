#!/usr/bin/env python

try:
    from setuptools import setup, find_packages
except ImportError:
    from distutils.core import setup
# from urllib import urlretrieve
# import tarfile, zipfile


setup(name='pevolution',
      version='1.0',
      description='This is a pipeline to find, align, and find trees for putatively related proteins.',
      author='Mihir Sarwade',
      author_email='mihir.sarwade@effem.com',
      packages=find_packages())

# urlretrieve("https://prank-msa.googlecode.com/files/prank.linux64.140110.tgz",'/tmp/prank.tgz')
# urlretrieve("http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip",'/tmp/PhyML.zip')
# urlretrieve("http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz",'/tmp/clustalw.tar.gz')
# tarfile.open('/tmp/prank.tgz').extractall(path='/usr/local/bin/')
# zipfile.ZipFile('/tmp/PhyML.zip').extractall(path='/usr/local/bin/')
# tarfile.open('/tmp/clustalw.tar.gz').extractall(path='.')
# tarfile.open('prottest-3.4-20140123.tar.gz').extractall(path='/usr/local/bin/')
