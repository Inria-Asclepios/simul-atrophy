#!/usr/bin/env python
import subprocess
import sys
import optparse

#ASSUMPTIONS:
#resultsPath must already exist and should contain the velocity file: vel.mha or output_stationary_velocity_field.mha
#patiendID should be present with the image to be warped.

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-i', '--inputFile', dest='inputFile', help='Filename of the input image to be warped')
parser.add_option('-d', '--displacementFile', dest='displacementFile', help='Filename of the field to be used to warp the input image')
parser.add_option('-o', '--outputFile', dest='outputFile', help='Filename of the output file')
parser.add_option('-v', '--invertField', dest='invertField', help='true or false: invert the field before warping or not')
parser.add_option('-n', '--number', dest='numOfScales', help='Total Number of time points to warp image to', type=int)
parser.add_option('-f', '--firstScale', dest='firstScale', help = 'Starting time point', type=int)
parser.add_option('-s', '--step', dest='stepSize', help = 'Step size between two time points', type=float)
parser.add_option('-c', '--compScheme', dest='compScheme', help = 'Computational scheme: 0 for scaling and squaring, is faster; 1 for forward euler, is more accurate')
(options, args) = parser.parse_args()

if options.inputFile is None:
    options.inputFile = raw_input('Enter the input image FileName')

print "your input File to be warped: " + options.inputFile

if options.displacementFile is None:
    options.displacementFile = raw_input('Enter the input field FileName')

if options.outputFile is None:
    options.outputFile = raw_input('Enter the output image FileName')

if options.invertField is None:
    options.invertField = raw_input('Enter invertField option: true or false')
    
if options.numOfScales is None:
    options.numOfScales = int(raw_input('Enter total number of time points'))

print "number of timePoints: " + str(options.numOfScales)

if options.firstScale is None:
    options.firstScale = int(raw_input('Enter starting time point:'))

print "first scale value: " + str(options.firstScale)

if options.stepSize is None:
    options.stepSize = float(raw_input('Enter step size:'))

print "step size: " + str(options.stepSize)

if options.compScheme is None:
    options.compScheme = (raw_input('Enter Computational scheme: 0 for scaling and squaring, is faster; 1 for forward euler, is more accurate'))

#executables
workDir = "/user/bkhanal/home/works/"
svfExp = workDir + "tools/svfExponential/SVFExponential"
warpImage = workDir + "AdLemModel/build/src/WarpImage"
combineImages = workDir + "tools/combineImages/build/combineImages"

#inputs of the executables
inputImage = options.inputFile
inputVel = options.displacementFile
scale = options.firstScale
num = 1
while num <= options.numOfScales:
    #outputs/results of the executables
    outputDisp = "dispScaledBy" + str(scale) + ".mha"
    outputImage = "warpedByDisp" + str(num) + ".mha"

    #create displacement field from inputVel
    svfExpOptions = " -s " + str(scale) + " -e " + options.compScheme + " -i " + inputVel + " -o " + outputDisp
    print svfExp + svfExpOptions + "\n"
    subprocess.call(svfExp + svfExpOptions, shell=True)

    #create warped Image from inputImage using outputDisp created with svfExp
    warpImageOptions = " --inImage " + inputImage + " --displacementImage " + outputDisp + " --outImage " + outputImage + " --invert " + options.invertField
    print warpImage + warpImageOptions + "\n"
    subprocess.call(warpImage + warpImageOptions, shell=True)

    #Delete displacement field
    print "rm " + outputDisp + "\n"
    subprocess.call("rm " + outputDisp, shell=True)

    #Go to next file
    num += 1
    scale += options.stepSize

#combine images: options: numOfScales imageType imageBaseFileName outputImageFileName
combineImagesOptions = " " + str(options.numOfScales) + " " + ".mha" + " " + "warpedByDisp" + " " + options.outputFile
print combineImages + combineImagesOptions + "\n"
subprocess.call(combineImages + combineImagesOptions, shell = True)
#delete inidvidual images:
#print "rm " + resultsPath + "warpedByDisp*.mha" + "\n"
#subprocess.call("rm " + resultsPath + "warpedByDisp*.mha", shell=True)

