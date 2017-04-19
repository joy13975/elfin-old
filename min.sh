#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Usage: min.sh <input_pdb> <max_cycles=200>"
	exit
fi

input="$1"
maxCycles="$2"
local="$3"

outDir=`dirname $input`
scOutput="${input/\.pdb/_min.sc}"

maxCycles=${maxCycles:-200}
local=${local:-"no"}
variant=${variant:-"default"}
release=${release:-"linuxgccrelease"}
wrapper=${wrapper:-""}

cmd="minimize.$variant.$release -overwrite -s $input -out:path:score $outDir -out:file:scorefile $scOutput -out:path:pdb $outDir -default_max_cycles $maxCycles"

if [[ "$local" == "yes" ]]; then
	$cmd
else
	sbatch -A other -p cpu -N 1 --ntasks-per-node=1 --wrap="$cmd"
fi
