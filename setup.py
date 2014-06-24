#!/usr/bin/env python

try:
  from setuptools import setup
except ImportError:
  from distutils.core import setup


setup(name='pevolution',
      version='1.0',
      description='This is a pipeline to find, align, and find trees for putatively related proteins.',
      author='Mihir Sarwade',
      author_email='mihir.sarwade@effem.com',
#      packages=['pevolution'],
      install_requires=['psutil >= 1.2.1','biopython >= 1.64'])

