#!/usr/bin/env python
import os
from sys import platform
from distutils.core import setup
from urllib.request import urlretrieve
import tarfile, zipfile, shutil


# setup(name='pevolution',
#       version='1.1',
#       description='This is a pipeline to find, align, and find trees for putatively related proteins.',
#       author='Mihir Sarwade',
#       author_email='mihir.sarwade@gmail.com',
#       #      packages=['pevolution'],
#       install_requires=['psutil >= 5.6.3','biopython >= 1.74', 'Bio'])

os.makedirs('tmp', exist_ok=True)
if platform == 'darwin':
    urlretrieve('http://wasabiapp.org/download/prank/prank.osx64.170427.zip','tmp/prank.zip')
    urlretrieve('http://phylipweb.github.io/phylip/download/phylip-3.695-osx.dmg','tmp/phylip.dmg')

    zipfile.ZipFile('tmp/prank.zip').extractall(path='tmp/')
    shutil.move('tmp/prank-msa', '.')

elif platform in ('win32', 'cygwin','win64'):
    urlretrieve('http://wasabiapp.org/download/prank/prank.cygwin.170427.zip',f'tmp{os.sep}prank.zip')
    urlretrieve('http://phylipweb.github.io/phylip/download/phylip-3.698.zip', f'tmp{os.sep}phylip.zip')
    urlretrieve('http://www.clustal.org/omega/clustal-omega-1.2.2-win64.zip', f'tmp{os.sep}clustal-omega.zip')
    zipfile.ZipFile(f'tmp{os.sep}prank.zip').extractall(path=f'tmp{os.sep}')
    # shutil.move(f'tmp{os.sep}prank-msa', '.')

else:
    urlretrieve('http://wasabiapp.org/download/prank/prank.linux.170427.tgz','tmp/prank.tgz')
    urlretrieve('http://wasabiapp.org/download/prank/prank.source.170427.tgz', 'tmp/prank.source.tgz')
    urlretrieve('http://phylipweb.github.io/phylip/download/phylip-3.697.tar.gz', 'tmp/phylip.tar.gz')

    tarfile.open('tmp/prank.tgz').extractall(path='tmp/prank-msa/')
    tarfile.open('tmp/phylip.tar.gz').extractall(path='tmp/phylip/')
    os.system('cd tmp/phylip/ && make -f Makefile.unx install')


urlretrieve("http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip", f'tmp{os.sep}PhyML.zip')
urlretrieve('https://github.com/ddarriba/prottest3/releases/download/3.4.2-release/prottest-3.4.2-20160508.tar.gz',
                f'tmp{os.sep}prottest.tar.gz')


zipfile.ZipFile(f'tmp{os.sep}PhyML.zip').extractall(path=f'tmp{os.sep}')
tarfile.open(f'tmp{os.sep}prottest.tar.gz').extractall(path=f'tmp{os.sep}')
