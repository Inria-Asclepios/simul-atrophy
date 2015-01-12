#!/usr/bin/env python
import subprocess
import sys
import optparse

workDir = "/epi/asclepios2/bkhanal/works/AdLemModel/"
target = workDir + "build/src/AdLemMain"

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-p', '--patientID', dest='patient', help='patient ID name which is also a valid directory containing the div/atrophy map')
parser.add_option('-a', '--atrophyFile', dest='atrophyFile',help='valid atrophy filename with extension present in patientID directory')
parser.add_option('-m', '--maskFile', dest='maskFile', help='valid brain Mask/segmentation filename with extension present in patientID data directory')
parser.add_option('--useTensorLambda', dest='useTensorLambda', action = 'store_true', help='If the option is provided, a valid lambdaFile must be given with the option -l')
parser.add_option('-l', '--lambdaFile', dest='lambdaFile', help='valid tensor image filename with extension present in patientID directory.')
parser.add_option('-f', '--resultsFilenamesPrefix', dest='resultsFilenamesPrefix', help='prefix to be added to all the output filenames')
parser.add_option('-j', '--jobFile', dest='jobFile', help='filename of the job pbs script.')
parser.add_option('-i', '--imageFile', dest='imageFile', help='valid input image that will be warped by the obtained displacement fields from the model.')
parser.add_option('-t', '--numOfTimeSteps', dest='numOfTimeSteps', help='(integer) number of time steps you want to solve the system.')
parser.add_option('--writePressure', dest='writePressure', action = 'store_true', help='If the option is not provided pressure map will not be written to the disk.')
parser.add_option('--writeForce', dest='writeForce', action = 'store_true', help='If the option is not provided, force term will not be written to the disk.')
parser.add_option('--writeResidual',dest='writeResidual', action = 'store_true', help='If the option is not provided, residual will not be written to the disk.')


(options, args) = parser.parse_args()
if options.patient is None:
    options.patient = raw_input('Enter a valid existing patient directory name')
resultsDir = workDir + "results/patients/" + options.patient + "/"

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
maskFile = resultsDir +  options.maskFile

if options.useTensorLambda is True:
    useTensorLambdaString = 'true'
    if options.lambdaFile is None:
        options.lambdaFile = raw_input('Enter a valid tensor image filename, or simply put a random string if you are not using tensor image by internally changing the PetscAdLemMain.cxx file')
    lambdaFile = resultsDir + options.lambdaFile
else:
    useTensorLambdaString = 'false'
    lambdaFile = "dummyFileName"

if options.resultsFilenamesPrefix is None:
    options.resultsFilenamesPrefix = raw_input('Enter a prefix string that will be added to all the output filenames')

if options.writePressure is True:
    writePressureString = 'true'
else:
    writePressureString = 'false'

if options.writeForce is True:
    writeForceString = 'true'
else:
    writeForceString = 'false'

if options.writeResidual is True:
    writeResidualString = 'true'
else:
    writeResidualString = 'false'

if options.jobFile is None:
    options.jobFile = raw_input('Enter a pbs script you want to call')

#Petsc options:
# These are here for reference purpose only since all these options will be explicitly provided in the .pbs script itself.
petscOptions = " -pc_type fieldsplit -pc_fieldsplit_type schur -pc_fieldsplit_schur_precondition self -pc_fieldsplit_dm_splits 0 -pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 "
#pc types for the fields
petscOptions += "-fieldsplit_0_pc_type hypre "
#monitor options
#petscOptions += "-fieldsplit_0_ksp_converged_reason -fieldsplit_0_ksp_max_it 100 -ksp_converged_reason"
petscOptions += "-fieldsplit_1_ksp_converged_reason -ksp_converged_reason "
#petscOptions += "-fieldsplit_1_ksp_max_it 3 -ksp_max_it 3 -ksp_rtol 1.0e-8 -ksp_monitor_true_residual"
petscOptions += "-fieldsplit_1_ksp_monitor_true_residual -ksp_monitor_true_residual -log_summary -ksp_view "

# petscCommand is also for reference purpose here, since will be explicitly provided in the .pbs script.
# petscCommand = "/opt/openmpi-gcc/current/bin/mpiexec "
# petscCommand += target + " -atrophyFile " + atrophyFile + " -useTensorLambda " + options.useTensorLambda + " -lambdaFile " + lambdaFile + " -maskFile " + maskFile + " -imageFile " + imageFile + " -numOfTimeSteps " + numOfTimeSteps + " -resPath " + resultsDir + " -resultsFilenamesPrefix " + options.resultsFilenamesPrefix + petscOptions

jobName = "AdLemModel" + options.patient
qsubCommand = "qsub " + workDir + "scripts/" + options.jobFile
qsubCommand += " -v atrophyFile=" + atrophyFile + ",useTensorLambda=" + useTensorLambdaString + ",lambdaFile=" + lambdaFile + ",maskFile=" + maskFile + ",imageFile=" + imageFile + ",numOfTimeSteps=" + options.numOfTimeSteps + ",resultsDir=" + resultsDir + ",resultsFilenamesPrefix=" + options.resultsFilenamesPrefix + ",writePressure=" + writePressureString + ",writeForce=" + writeForceString + ",writeResidual=" + writeResidualString +  " -d " + resultsDir + " -N " + jobName
print "given qsub command: \n" + qsubCommand
subprocess.call(qsubCommand, shell=True)




