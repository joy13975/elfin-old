#!/bin/bash

bmDir=${bmDir:-"l10"}
outDir=$bmDir"_pdb"

mkdir -p bm/$outDir
for dir in `ls bm/$bmDir`; do 
	if [[ -d "bm/$bmDir/$dir" ]]; then 
		for f in `ls bm/$bmDir/$dir/*.pdb`; do 
			cp $f bm/$outDir/$dir".pdb" 
			break 1			
		done; 
	fi; 
done
