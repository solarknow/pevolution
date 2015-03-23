pevolution
==========

This is a pipeline to find, align, and find trees for putatively related proteins.

----------
Before installing make sure you have your email set as ENTREZ_EMAIL as an environmental variable.
To install pre-reqs,
--------------------
##Debian/Ubuntu:

```shell

sudo apt-get install python-pip mrbayes mrbayes-doc wget python-dev unzip ncbi-blast+ clustalw openjdk-7-jre paml prank

tar zxvf prottest-3.4-20140123.tar.gz 
mv prottest-3.4-20140123 prottest

python setup.py install
```
------------------
##RHEL/Fedora/CentOS:

Run ./installdeps.rhel
as root.

--------
##Mac OS X:

If you don't have python installed,
Install python: https://www.python.org/ftp/python/2.7.6/python-2.7.6-macosx10.6.
dmg
Install pip: https://pip.pypa.io/en/latest/installing.html
Install MrBayes: http://downloads.sourceforge.net/project/mrbayes/mrbayes/3.2.2/MrBayes-3.2.2_installer_MACx64.pkg

To install dependencies, run sudo python setup.py install in this directory.
-------
##Windows:

If you don't have python installed,
Install python: https://www.python.org/ftp/python/2.7.6/python-2.7.6.msi
Make sure the directory you installed Python in (usually C:\Python27) is in the 
PATH. 
Install pip: https://pip.pypa.io/en/latest/installing.html according to the python you have installed (in the above step, python 2.7 32bit)

Open a cmd.exe in this directory and run: python setup.py install
