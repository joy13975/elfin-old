#!/bin/bash

for dir in bm/funOut/*; do
	for f in $dir/*; do
		echo $f
		break 1
	done
done
