#!/bin/bash

bmDir=${bmDir:-"l10"}
outDir=$bmDir"_nat"

mkdir -p bm/$outDir
for dir in `ls bm/$bmDir`; do 
	if [[ -d "bm/$bmDir/$dir" ]]; then 
		for f in `ls bm/$bmDir/$dir/*.pdb`; do 
			if [[ $f == *"_0001.pdb" ]]; then
				continue
			fi

			cp $f bm/$outDir/$dir".pdb" 
			break 1			
		done; 
	fi; 
done
