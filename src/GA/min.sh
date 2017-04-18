#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: min.sh <input_pdb> <max_cycles=200>"
	exit
fi

input="$1"
outDir=`dirname $input`
scOutput="${input/\.pdb/_min.sc}"
maxCycles="$2"

maxCycles=${maxCycles:-200}

cmd="minimize.default.linuxgccrelease -overwrite -s $input -out:path:score $outDir -out:file:scorefile $scOutput -out:path:pdb $outDir -default_max_cycles $maxCycles"
#$cmd
sbatch -A other -p cpu -N 1 --ntasks-per-node=1 --wrap="$cmd"
