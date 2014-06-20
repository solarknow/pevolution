pevolution
==========

This is a pipeline to find, align, and find trees for putatively related proteins.

----------
To install pre-reqs, 
Debian/Ubuntu:
sudo apt-get install python-psutil python-biopython mrbayes mrbayes-doc 
wget https://prank-msa.googlecode.com/files/prank.linux64.140110.tgz
tar zxvf prank.linux64.140110.tgz prank
wget http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip
unzip PhyML-3.1.zip phyml

RHEL/Fedora/CentOS:
sudo yum install python-psutil python-biopython libbeagle-dev
wget https://prank-msa.googlecode.com/files/prank.linux64.140110.tgz
tar zxvf prank.linux64.140110.tgz prank
wget http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip
unzip PhyML-3.1.zip phyml
wget http://downloads.sourceforge.net/project/mrbayes/mrbayes/3.2.2/mrbayes-3.2.2.tar.gz
tar zxvf mrbayes-3.2.2.tar.gz
cd mrbayes_3.2.2


Mac OS X:
If you don't have python installed,
Install python: https://www.python.org/ftp/python/2.7.6/python-2.7.6-macosx10.6.
dmg
Install pip: https://pip.pypa.io/en/latest/installing.html
Install MrBayes: http://downloads.sourceforge.net/project/mrbayes/mrbayes/3.2.2/MrBayes-3.2.2_installer_MACx64.pkg

To install dependencies, run sudo python setup.py install in this directory
ory.

Windows:
If you don't have python installed,
Install python: https://www.python.org/ftp/python/2.7.6/python-2.7.6.msi
Make sure the directory you installed Python in (usually C:\Python27) is in the 
PATH. 
Install pip: https://pip.pypa.io/en/latest/installing.html according to the python you have installed (in the above step, python 2.7 32bit)

Open a cmd.exe in this directory and run: python setup.py install
