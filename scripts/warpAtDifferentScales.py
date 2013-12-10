#!/usr/bin/env python
import subprocess
import sys
import optparse

#ASSUMPTIONS:
#resultsPath must already exist and should contain the velocity file: vel.mha or output_stationary_velocity_field.mha
#patiendID should be present with the image to be warped.

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-p', '--patientID', dest='patientId', help='ID of the patient, a valide directory name that contains the source brain image to be warped.')
parser.add_option('-i', '--inputFile', dest='inputFile', help='Filename of the input image to be warped')
parser.add_option('-n', '--number', dest='numOfScales', help='Total Number of time points to warp image to', type=int)
parser.add_option('-f', '--firstScale', dest='firstScale', help = 'Starting time point', type=int)
parser.add_option('-s', '--step', dest='stepSize', help = 'Step size between two time points', type=float)
parser.add_option('-c', '--choice', dest='choice', help = 'r for using velocity computed from registration; b for using velocity computed from biophyisical model.')

(options, args) = parser.parse_args()
if options.patientId is None:
    options.patientId = raw_input('Enter patient ID, a directory that contains the source image to be warped')

print "your patient: " + options.patientId

if options.inputFile is None:
    options.inputFile = raw_input('Enter the input image FileName')

print "your input File to be warped: " + options.inputFile

if options.numOfScales is None:
    options.numOfScales = int(raw_input('Enter total number of time points'))

print "number of timePoints: " + str(options.numOfScales)

if options.firstScale is None:
    options.firstScale = int(raw_input('Enter starting time point:'))

print "first scale value: " + str(options.firstScale)

if options.stepSize is None:
    options.stepSize = float(raw_input('Enter step size:'))

print "step size: " + str(options.stepSize)

print "registration vs. biophysical model choice: " + options.choice

while ( not ((options.choice == "r") or (options.choice == "b")) ):
    options.choice = raw_input('Enter valid choice for velocity field to use: r: from registration; b: from biophysical model')

#base paths
workDir = "/user/bkhanal/home/works/"

if options.choice == "r":
    resultsPath = "/user/bkhanal/home/data/patients/" + options.patientId + "/"
    velName = "output_stationary_velocity_field.mha"
else:
    resultsPath = workDir + "AdLemModel/results/patients/" + options.patientId + "/"
    velName = "vel.mha"

#executables
svfExp = workDir + "tools/svfExponential/SVFExponential"
warpImage = workDir + "tools/warpImage/build/warpImage"
combineImages = workDir + "tools/combineImages/build/combineImages"

#inputs of the executables
inputImage = "/user/bkhanal/home/data/patients/" + options.patientId + "/" + options.inputFile
inputVel = resultsPath + velName

scale = options.firstScale

num = 1
while num <= options.numOfScales:
    #outputs/results of the executables
    outputDisp = resultsPath + "dispScaledBy" + str(scale) + ".mha"
    outputImage = resultsPath + "warpedByDisp" + str(num) + ".mha"

    #create displacement field from inputVel
    svfExpOptions = " -s " + str(scale) + " -i " + inputVel + " -o " + outputDisp
    print svfExp + svfExpOptions + "\n"
    subprocess.call(svfExp + svfExpOptions, shell=True)

    #create warped Image from inputImage using outputDisp created with svfExp
    warpImageOptions = " " + inputImage + " " + outputDisp + " " + outputImage
    print warpImage + warpImageOptions + "\n"
    subprocess.call(warpImage + warpImageOptions, shell=True)

    #Delete displacement field
    print "rm " + outputDisp + "\n"
    subprocess.call("rm " + outputDisp, shell=True)

    #Go to next file
    num += 1
    scale += options.stepSize

#combine images: options: numOfScales imageType imageBaseFileName outputImageFileName
if options.choice == "r":
    outputImageFileName = "image4dFromRegistration.mha"
else:
    outputImageFileName = "image4dFromBiophysicalModel.mha"

combineImagesOptions = " " + str(options.numOfScales) + " " + ".mha" + " " + resultsPath + "warpedByDisp" + " " + resultsPath + outputImageFileName
print combineImages + combineImagesOptions + "\n"
subprocess.call(combineImages + combineImagesOptions, shell = True)
#delete inidvidual images:
print "rm " + resultsPath + "warpedByDisp*.mha" + "\n"
subprocess.call("rm " + resultsPath + "warpedByDisp*.mha", shell=True)

