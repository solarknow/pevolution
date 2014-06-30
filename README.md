pevolution
==========

This is a pipeline to find, align, and find trees for putatively related proteins.

----------
Before installing make sure you have your email set as ENTREZ_EMAIL as an environmental variable.
To install pre-reqs,
--------------------
##Debian/Ubuntu:

```shell

sudo apt-get install python-pip mrbayes wget python-dev unzip ncbi-blast+ openjdk-7-jre prank phyml

sudo python setup.py install
```
------------------
##RHEL/Fedora/CentOS:

```shell

sudo yum install python-pip python-devel java-1.7.0-openjdk-devel wget autoconf automake svn unzip ncbi-blast+
sudo yum groupinstall "Development Tools"
```
###Install libhmsbeagle for mrbayes

```shell

svn checkout http://beagle-lib.googlecode.com/svn/trunk/ beagle-lib
cd beagle-lib
./autogen.sh
./configure
make
sudo make install
```
###Install MrBayes

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
###Other Deps (no compilation necessary)

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
