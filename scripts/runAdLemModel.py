#!/usr/bin/env python
import subprocess
import sys
import optparse

workDir = "/user/bkhanal/home/works/AdLemModel/"
target = workDir + "build/src/AdLemMain"
petscCommand = "/home/bkhanal/Documents/softwares/petsc-3.4.3/bin/petscmpiexec -n 4 "

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-p', '--patientID', dest='patient', help='patient ID name which is also a valid directory containing the div/atrophy map')
parser.add_option('-a', '--atrophyFile', dest='atrophyFile',help='valid atrophy filename with extension present in patientID directory')
parser.add_option('-m', '--maskFile', dest='maskFile', help='valid brain Mask/segmentation filename with extension present in patientID data directory')
parser.add_option('-u', '--useTensorLambda', dest='useTensorLambda', help='provide true or false. If true given, must provide a valid lambdaFile with the option -l')
parser.add_option('-l', '--lambdaFile', dest='lambdaFile', help='valid tensor image filename (if useTensorLambda is set to true) with extension present in patientID directory.')
parser.add_option('-f', '--resultsFilenamesPrefix', dest='resultsFilenamesPrefix', help='prefix to be added to all the output filenames')
parser.add_option('-i', '--imageFile', dest='imageFile', help='valid input image that will be warped by the obtained displacement fields from the model.')
parser.add_option('-t', '--numOfTimeSteps', dest='numOfTimeSteps', help='(integer) number of time steps you want to solve the system.')

(options, args) = parser.parse_args()

if options.patient is None:
    options.patient = raw_input('Enter a valid existing patient directory name')
#Data and Results Directory:
resultsDir = workDir + "results/patients/" + options.patient + "/"
#create if not created: Nah, should be already created to have atrophyFile and MaskFile.
#subprocess.call("mkdir -p " + resultsDir, shell=True) 

if options.numOfTimeSteps is None:
    options.numOfTimeSteps = raw_input('Enter a number of time steps you want to solve the system.')

if options.imageFile is None:
    options.imageFile = raw_input('Enter a valid existing image file that will be warped by the obtained displacement fields')
imageFile = resultsDir + options.imageFile

if options.atrophyFile is None:
    options.atrophyFile = raw_input('Enter a valid atrophy filename')
atrophyFile = resultsDir + options.atrophyFile

if options.maskFile is None:
    options.maskFile = raw_input('Enter a valid mask/segmented image filename')
maskFile = resultsDir + options.maskFile

if options.useTensorLambda is None:
    options.useTensorLambda = raw_input('provide true or false depending on if you intend to use a lambda from a tensor image.')

if options.useTensorLambda is True:
    if options.lambdaFile is None:
        options.lambdaFile = raw_input('Enter a valid tensor image filename, or simply put a random string if you are not using tensor image by internally changing the PetscAdLemMain.cxx file')
    lambdaFile = resultsDir + options.lambdaFile
else:
    lambdaFile = "dummyFileName"

if options.resultsFilenamesPrefix is None:
    options.resultsFilenamesPrefix = raw_input('Enter a prefix string that will be added to all the output filenames')


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
petscOptions += "-fieldsplit_1_ksp_monitor_true_residual -ksp_monitor_true_residual" # -log_summary" # -ksp_view"


petscCommand += target + " -atrophyFile " + atrophyFile + " -useTensorLambda " + options.useTensorLambda + " -lambdaFile " + lambdaFile + " -maskFile " + maskFile + " -imageFile " + imageFile + " -numOfTimeSteps " + options.numOfTimeSteps + " -resPath " + resultsDir + " -resultsFilenamesPrefix " + options.resultsFilenamesPrefix + petscOptions

print "given command:\n" + petscCommand + "\n"
subprocess.call(petscCommand, shell=True)




