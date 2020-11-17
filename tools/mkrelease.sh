#!/bin/sh

ver=2.0.0
workdir=./minfft-$ver

rm -rf $workdir $workdir.tar.gz
mkdir $workdir
cp ../README ../LICENSE ../minfft.h ../minfft.c ../minfft.F03 $workdir
tar cf $workdir.tar $workdir --owner=0 --group=0
gzip -9 $workdir.tar
echo Created $workdir.tar.gz
