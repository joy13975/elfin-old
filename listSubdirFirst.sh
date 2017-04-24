#!/bin/bash

dir=${dir:-"bm/l10/"}

extension=${extension:-".json"}

for subdir in $dir/*; do
	if [[ -d "$subdir" ]]; then
		for f in $subdir/*$extension; do
			echo $f
			break 1
		done
	fi
done