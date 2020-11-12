#!/bin/bash

mkdir problems
cd problems

# was missing.. idk
cpanm  Hash::StoredIterator

export PERL_MM_USE_DEFAULT=1
cpan Array::Utils

# for solving IntelliJ IDEA problem with connecting to debugger ...
git clone https://github.com/Camelcade/Devel-Camelcadedb
cd Devel-Camelcadedb

# TODO in 'lib/Devel/Camelcadedb.pm' change 'our $VERSION = "v<CHANGE>"' 

perl Makefile.PL
make
make test
make install

