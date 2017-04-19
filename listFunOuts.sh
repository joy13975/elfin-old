#!/bin/bash

extension=${extension:-".json"}

for dir in bm/funOut/*; do
	for f in $dir/*$extension; do
		echo $f
		break 1
	done
done
