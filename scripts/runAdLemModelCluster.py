#!/usr/bin/env python
import subprocess
import sys
import optparse

workDir = "/epi/asclepios2/bkhanal/works/AdLemModel/"
dataDir = "/epi/asclepios2/bkhanal/data/patients/"
target = workDir + "build/src/AdLemMain"


#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-p', '--patientID', dest='patient', help='patient ID name which is also a valid directory containing the div/atrophy map')

(options, args) = parser.parse_args()
if options.patient is None:
    options.patient = raw_input('Enter a valid existing patient directory name')

#Input Selection:
divFile = dataDir + options.patient + "/divergence.nii"
#results Directory:
resultsDir = workDir + "results/patients/" + options.patient + "/"
#create if not created:
subprocess.call("mkdir -p " + resultsDir, shell=True)

#Petsc options:
petscOptions = " -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_dm_splits 0 -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 "
#pc types for the fields
petscOptions += "-fieldsplit_0_pc_type hypre "
#monitor options
petscOptions += "-fieldsplit_0_ksp_converged_reason -fieldsplit_0_ksp_max_it 100 -ksp_converged_reason"


petscCommand = "/opt/openmpi-gcc/current/bin/mpiexec "
petscCommand += target + " -divFile " + divFile + " -resPath " + resultsDir + petscOptions

jobName = "AdLemModel" + options.patient
qsubCommand = "qsub " + workDir + "scripts/AdLem3Dhypre1.pbs -v divFile=" + divFile + ",resultsDir=" + resultsDir + " -d " + resultsDir + " -N " + jobName
print "given qsub command: \n" + qsubCommand
subprocess.call(qsubCommand, shell=True)




