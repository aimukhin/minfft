#!/bin/sh

ver=1.2.0
workdir=./minfft-$ver

rm -rf $workdir $workdir.tar.gz
mkdir $workdir
cp ../README ../LICENSE ../minfft.{h,c,F03} $workdir
tar cf $workdir.tar $workdir --owner=0 --group=0
gzip -9 $workdir.tar
echo Created $workdir.tar.gz
