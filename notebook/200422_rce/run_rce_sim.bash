#!/bin/bash

panfa=cov1.bylength.fa

if [ ! -f $panfa ] ; then
	echo Not found panfa=$panfa >> /dev/stderr
	exit 1
fi

NrReads=10000

for ReadLength in 100 150
do
	for PctId in 100 95 90 85 80
	do
		echo === $PctId L$ReadLength ===
		simple_paired_read_simulator.py $panfa \
		  $ReadLength \
		  $NrReads \
		  $PctId \
		  sim${PctId}_${ReadLength}_1.fq \
		  sim${PctId}_${ReadLength}_2.fq
	done
done
