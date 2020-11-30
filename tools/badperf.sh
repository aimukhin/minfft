#!/bin/sh

if [ -n "$1" ]; then
	t=$1
else
	t=10
fi

while read x s a p1 p2; do
	r=$((100*(p1-p2)/p2))
	if [ $r -lt -$t ]; then
		echo "$r	$x $s"
	fi
done
