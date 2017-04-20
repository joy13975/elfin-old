#!/bin/bash

trap "exit" INT

if [[ $# < 1 ]]; then
	echo "Usage: CloseLoopAll.sh <input>"
	exit
fi

input=$1
input_dir=`dirname $input`
variant=${variant:-"default"}
release=${release:-"linuxgccrelease"}
local=${local:="no"}

cmd="loopmodel.$variant.$release -in:file:s $input -loops:loop_file ${input/\.pdb/.loops} -out:path:all $input_dir @res/loopModel.config"
# echo CMD is $cmd

if [[ "$local" == "yes" ]]; then
	$cmd
else
	sbatch -A other -p cpu -N 1 --ntasks-per-node=1 --wrap="$cmd"		
fi

