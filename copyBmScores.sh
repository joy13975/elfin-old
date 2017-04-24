#!/bin/bash

bmDir=${bmDir:-"l10"}
outDir=$bmDir"_sc"

mkdir -p bm/$outDir
for dir in `ls bm/$bmDir`; do 
	if [[ -d "bm/$bmDir/$dir" ]]; then 
		for f in `ls bm/$bmDir/$dir/*_comp.sc`; do 
			cp $f bm/$outDir/$dir"_comp.sc" 
			break 1			
		done; 
	fi; 
done
