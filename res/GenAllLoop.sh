#!/bin/bash

trap "exit" INT

if [[ $# < 1 ]]; then
	echo 'Usage: ./GenAllLoop.sh <input_dir>'
	exit
fi

input_dir=$1

for file in $input_dir/*_mc_*.pdb; do
	../src/Python/GenLoopRegion.py --input $file > ${file/\.pdb/.loops}
done
