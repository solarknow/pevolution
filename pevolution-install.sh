#!/bin/bash

apt-get update
apt-get install -y python-pip wget python-dev unzip default-jre subversion
cpan install LWP:Simple

apt-get install -y ncbi-blast+ clustalw phyml mrbayes mrbayes-doc

python setup.py install

tar zxvf prottest-3.4.2-20160508.tar.gz
