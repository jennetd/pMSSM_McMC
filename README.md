#Markov chain Monte Carlo code for scanning the 19D pMSSM parameter space
Original code from Malte Mrowietz

## Docker image 
An image containing the necessary packages is available: jdickins/pmssm-env:latest

To run the docker container on your personal computer:
'''
docker pull jdickins/pmssm-env:latest
docker run -it jdickins/pmssm-env:latest /bin/bash
'''

To run the docker container using singularity on lxplus:
'''
singularity exec -B /afs /cvmfs/unpacked.cern.ch/registry.hub.docker.com/jdickins/pmssm-env:latest /bin/bash
'''

To run a test scan, run the following inside the container:
'''
source /root-6.12.06-build/bin/thisroot.sh
git clone https://github.com/jennetd/pMSSM_McMC workdir
cd pMSSM_McMC
'''

External (public) packages are contained in the packages/ directory:
 * SPheno 4.0.4 - can use pre-compiled version in docker container
 * FeynHiggs 2.16.1 - can use pre-compiled version in docker container
 * superiso 4.0 - can use pre-compiled version in docker container
 * GM2Calc 1.7.3 - can use pre-compiled version in docker container
 * HiggsBounds 5.9.1 - need to install locally (because the contents must be editable by scan job)
 * HiggsSignals 2.6.0 - need to install locally (because the contents must be editable by scan job)
 * Micromegas 5.2.4 - not needed at this time

## Running the code

To install HiggsSignals and HiggsBounds, do
'''
./install-hs.sh
./install-hb.sh
'''

To run an example scan over 50 points:
'''
mkdir output
python mcmc.py -m new -n 50 -o ouptut
'''

The full list of options :
* `-m, --mode ["new","resume"]`:  
   new: starts a new chain from a random point in the pMSSM, sampling from a flat prior in the 19 pMSSM dimensions.  
    resume: resumes a previously created chain. Specify the root file containing the previous chain. The chain is continued from the point in the input root file that has the highest iteration index. The input file is specified with the `-i, [input]` argument.
* `-i, --input [input]`,If runmode resume, specify input root file (default = None)
* `-o, --output [output]`,Specify an output directory (required)
* `-n, --npoints [npoints]`, How many points to run the MCMC for (default=1000)
* `-c, --chain_index [index]`, chain index for the chain. Points are later identified by the chain_index and the iteration index (default = 1)
* `-mi, --move_interval [interval]` , How many points to generate before starting a new root file and moving the already created points to the output directory (default = 300)

## Editing the code

Adding observables:
 * In likelihood.py: add the value for the observable and its uncertainty to the `likelihood_contributions` dictionary. To change the way the likelihood is computed for the McMC, edit likelihood.py
 * In mcmc.py: add a way to retrieve the observable to the `get_observables` function. This may require you to write an interface to external tools.  

Adjusting parameter ranges:
 * In mcmc.py: edit the tuple `parameter_ranges["parameter"]` at the top of the file

Adjusting scan parameters:
 * In mcmc.py: the base of the logarithm used for log stepping is stored in the variable `base`. The width of the gaussian that determines stepping is: `log(parameter_max - parameter_min) * width_coefficient`. Change this by changing the value of `width_coefficient`. 

## Running batch jobs

To run HTCondor jobs on lxplus, visit the subdirectory lxplus-condor. Create a tarball of the code version you want to run. Update the name of the tarball in submit-scan.sh and run-scan.sh. Edit the job name and number of points to scan in submit-scan.sh.  Make sure the correct tarball is copied to the remote working directory by run_scan.sh.

Submit the job by running
'''
./submit-scan.sh
'''