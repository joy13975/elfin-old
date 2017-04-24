#!/bin/bash

bmDir=${bmDir:-"l10"}
outDir=$bmDir"_relaxed"

mkdir -p bm/$outDir
for dir in `ls bm/$bmDir`; do 
	if [[ -d "bm/$bmDir/$dir" ]]; then 
		for f in `ls bm/$bmDir/$dir/*_0001.pdb`; do 
			cp $f bm/$outDir/$dir"_0001.pdb" 
			break 1			
		done; 
	fi; 
done
