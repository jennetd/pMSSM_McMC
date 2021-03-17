#!/bin/bash

npoints=100
name="test9_"${npoints}

subfile=logs/${name}.sub
rm $subfile

echo "universe                = docker" >> $subfile
echo "docker_image            = jdickins/pmssm-env:latest" >> $subfile
echo "executable              = run-scan.sh" >> $subfile
echo "arguments               = ${npoints}" >> $subfile
echo "should_transfer_files   = YES" >> $subfile
echo "when_to_transfer_output = ON_EXIT" >> $subfile
echo "transfer_input_files    = pMSSM_McMC-master.tar.gz" >> $subfile
echo "transfer_output_files   = pMSSM_McMC/output_${npoints}" >> $subfile
echo "output                  = logs/${name}.out" >> $subfile
echo "error                   = logs/${name}.err" >> $subfile
echo "log                     = logs/${name}.log" >> $subfile

# Job flavour determines job wall time
# https://batchdocs.web.cern.ch/local/submit.html#job-flavours
echo "+JobFlavour             = \"longlunch\"" >> $subfile

echo "queue" >> $subfile

condor_submit $subfile
