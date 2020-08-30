pevolution
==========

This is a pipeline to find, align, and find trees for putatively related proteins.

----------

Before installing make sure you have your email set as ENTREZ_EMAIL as an environmental variable.
To install pre-reqs,

--------------------

## Debian/Ubuntu:

```shell

sudo apt-get install python-pip mrbayes mrbayes-doc wget python-dev unzip ncbi-blast+ clustalw openjdk-7-jre paml

tar zxvf prottest-3.4-20140123.tar.gz 
mv prottest-3.4-20140123 prottest

python setup.py install
```
------------------

## RHEL/Fedora/CentOS:

```shell

sudo yum install python-pip python-devel java-1.7.0-openjdk-devel wget autoconf automake svn unzip ncbi-blast+
sudo yum groupinstall "Development Tools"
```
### Install libhmsbeagle for mrbayes

```shell

svn checkout http://beagle-lib.googlecode.com/svn/trunk/ beagle-lib
cd beagle-lib
./autogen.sh
./configure
make
sudo make install
```
### Install MrBayes

```shell

wget http://downloads.sourceforge.net/project/mrbayes/mrbayes/3.2.2/mrbayes-3.2.2.tar.gz
tar zxvf mrbayes-3.2.2.tar.gz
cd mrbayes_3.2.2
autoconf
./configure
make
sudo make install
```
---------

### Other Deps (no compilation necessary)

```shell
wget https://prank-msa.googlecode.com/files/prank.linux64.140110.tgz
tar zxvf prank.linux64.140110.tgz prank
wget http://www.atgc-montpellier.fr/download/binaries/phyml/PhyML-3.1.zip
unzip PhyML-3.1.zip
mv Phyml-3.1 phyml
wget http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
tar zxvf clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
mv clustalw-2.1-linux-x86_64-libcppstatic clustalw
sudo python setup.py install
```
--------

## Mac OS X:

If you don't have python installed,
Install python: `brew install python@3.8`
Install blast: `brew install blast`
Install MrBayes: `wget https://github.com/NBISweden/MrBayes/releases/download/v3.2.7/mrbayes-3.2.7.tar.gz && tar zxvf mrbayes-3.2.7.tar.gz`
Install PRANK:  
```wget http://wasabiapp.org/download/prank/prank.osx64.170427.tgz && tar zxvf prank.osx64.170427.tgz```  
Install prottest:
`wget https://github.com/ddarriba/prottest3/releases/download/3.4.2-release/prottest-3.4.2-20160508.tar.gz && tar zxvf prottest-3.4.2-20160508.tar.gz`  
To install python dependencies, run sudo python3 setup.py install in this directory.

-------

## Windows:

If you don't have python installed,
Install python: https://www.python.org/ftp/python/2.7.6/python-2.7.6.msi
Make sure the directory you installed Python in (usually C:\Python27) is in the 
PATH. 
Install pip: https://pip.pypa.io/en/latest/installing.html according to the python you have installed (in the above step, python 2.7 32bit)

Open a cmd.exe in this directory and run: python setup.py install
