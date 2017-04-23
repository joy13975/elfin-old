#!/bin/bash

extension=${extension:-".json"}

for dir in bm/l10/*; do
	if [ -d $dir ]; then
		for f in $dir/*$extension; do
			echo $f
			break 1
		done
	fi
done
