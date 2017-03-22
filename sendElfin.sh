#!/bin/bash

: ${elfin}
: ${queueName="batch"}
: ${ppn="256:knl"}
: ${runScript="src/GA/run.sh"}
: ${configFile="src/GA/settings.json"}
: ${outputDir="output"}
: ${walltime="00:30:00"}


qsub -N $jobName \
	-q $queueName \
	-joe -o $outputDir/log \
	 -l nodes=1:ppn=$ppn,walltime=$walltime \
	 $runScript