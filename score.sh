#!/bin/bash

trap "exit" INT

if [ $# -lt 1 ]; then
	echo "Usage: score.sh <min_pdb> <native_pdb>"
	exit
fi

min="$1"
native="$2"

scOutput="${native/\.pdb/_comp.sc}"

local=${local:-"no"}
variant=${variant:-"default"}
release=${release:-"linuxgccrelease"}
wrapper=${wrapper:-""}

cmd="score.$variant.$release -overwrite -in:file:native $native -s $min -out:file:scorefile $scOutput"

if [[ "$local" == "yes" ]]; then
	$cmd
else
	sbatch -A other -p cpu -N 1 --ntasks-per-node=1 --wrap="$cmd"
fi
