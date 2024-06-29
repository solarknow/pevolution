pevolution
==========

This is a pipeline to find, align, and find trees for putatively related proteins.

----------

Before installing make sure you have your email set as ENTREZ_EMAIL as an environmental variable.
To install pre-reqs,

--------------------

## Debian/Ubuntu:

```shell

sudo apt-get install python3-pip mrbayes mrbayes-doc wget python3-dev ncbi-blast+ clustalw openjdk-8-jre paml prank

tar zxvf prottest-3.4.2-20160508.tar.gz 
mv prottest-3.4.2 prottest

python setup.py install
```
------------------

## RHEL/Fedora/CentOS:

```shell

sudo yum install python3-pip python3-devel java-1.8.0-openjdk-devel wget autoconf automake svn ncbi-blast+
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

git clone --depth=1 https://github.com/NBISweden/MrBayes.git
cd MrBayes
./configure
make
sudo make install
```
---------

### Other Deps (no compilation necessary)

```shell
wget http://wasabiapp.org/download/prank/prank.linux64.170427.tgz
tar zxvf prank.linux64.170427.tgz prank
wget https://github.com/stephaneguindon/phyml/archive/v3.3.20190321.tar.gz
tar zxvf v3.3.20190321.tar.gz
mv phyml-3.3.20190321 phyml
wget http://www.clustal.org/download/current/clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
tar zxvf clustalw-2.1-linux-x86_64-libcppstatic.tar.gz
mv clustalw-2.1-linux-x86_64-libcppstatic clustalw
sudo python setup.py install
```