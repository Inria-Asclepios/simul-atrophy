#!/usr/bin/env python
import subprocess
import sys
import optparse

workDir = "/user/bkhanal/home/works/AdLemModel/"
dataDir = "/user/bkhanal/home/data/patients/"
target = workDir + "build/src/AdLemMain"


#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-p', '--patientID', dest='patient', help='patient ID name which is also a valid directory containing the div/atrophy map')
parser.add_option('-a', '--atrophyFile', dest='divFile',help='valid atrophy/divergence filename with extension present in patientID directory')
parser.add_option('-m', '--maskFile', dest='maskFile', help='valid brain Mask/segmentation filename with extension present in patientID data directory')

(options, args) = parser.parse_args()
if options.patient is None:
    options.patient = raw_input('Enter a valid existing patient directory name')
if options.divFile is None:
    options.divFile = raw_input('Enter a valid div or atrophy filename')
if options.maskFile is None:
    options.maskFile = raw_input('Enter a valid mask/segmented image filename')

#Input Selection:
#divFile = dataDir + options.patient + "/divergence.nii"
divFile = dataDir + options.patient + "/" + options.divFile
maskFile = dataDir + options.patient + "/" + options.maskFile
#results Directory:
resultsDir = workDir + "results/patients/" + options.patient + "/"
#create if not created:
subprocess.call("mkdir -p " + resultsDir, shell=True)

#Petsc options:
petscOptions = " -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_precondition self -pc_fieldsplit_dm_splits 0 -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 "
#petscOptions += "-pc_fieldsplit_schur_fact_type upper "
#petscOptions = " -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_detect_saddle_point "
#pc types for the fields
petscOptions += "-fieldsplit_0_pc_type hypre "
#petscOptions += "-fieldsplit_0_pc_type lu -fieldsplit_0_pc_factor_mat_solver_package superlu_dist "
#monitor options
petscOptions += "-fieldsplit_1_ksp_converged_reason -ksp_converged_reason "
#petscOptions += "-fieldsplit_1_ksp_max_it 3 -ksp_max_it 3 "
petscOptions += "-fieldsplit_1_ksp_monitor_true_residual -ksp_monitor_true_residual" # -ksp_view"

petscCommand = "/home/bkhanal/Documents/softwares/petsc-3.4.3/bin/petscmpiexec -n 4 "
petscCommand += target + " -divFile " + divFile + " -maskFile " + maskFile + " -resPath " + resultsDir + petscOptions

print "given command:\n" + petscCommand + "\n"
subprocess.call(petscCommand, shell=True)




