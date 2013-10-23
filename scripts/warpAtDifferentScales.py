#!/usr/bin/env python
import subprocess
import sys
import optparse

#base paths
workDir = "/user/bkhanal/home/works/"
resultsPath = workDir + "AdLemModel/results/"
dataPath = workDir + "AdLemModel/data/"

#executables
svfExp = workDir + "tools/svfExponential/SVFExponential"
warpImage = workDir + "tools/warpImage/build/warpImage"
combineImages = workDir + "tools/combineImages/build/combineImages"


#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-i', '--inputFile', dest='inputFile', help='Filename of the input image to be warped')
parser.add_option('-v', '--velocityFile', dest='velFile', help='Filename of the input velocity used to create displacements')
parser.add_option('-n', '--number', dest='numOfScales', help='Total Number of time points to warp image to', type=int)
parser.add_option('-f', '--firstScale', dest='firstScale', help = 'Starting time point', type=int)
parser.add_option('-s', '--step', dest='stepSize', help = 'Step size between two time points', type=int)
parser.add_option('-o', '--outputFile', dest='outputFile', help='Filename of the output 4D image that is warped at sequence of times')

(options, args) = parser.parse_args()
if options.inputFile is None:
    options.inputFile = raw_input('Enter the input image FileName')

if options.velFile is None:
    options.velFile = raw_input('Enter the velocity image FileName')

if options.numOfScales is None:
    options.numOfScales = int(raw_input('Enter total number of time points'))

if options.firstScale is None:
    options.firstScale = int(raw_input('Enter starting time point:'))

if options.stepSize is None:
    options.stepSize = int(raw_input('Enter step size:'))

if options.outputFile is None:
    options.outputFile = raw_input('Enter the output image FileName')


#inputs of the executables
inputImage = dataPath + "Template/" + options.inputFile
inputVel = resultsPath + options.velFile

scale = options.firstScale

num = 1
while num <= options.numOfScales:
    #outputs/results of the executables
    outputDisp = resultsPath + "dispScaledBy" + str(scale) + ".mha"
    outputImage = resultsPath + "warpedByDisp" + str(num) + ".mha"

    #create displacement field
    svfExpOptions = " -s " + str(scale) + " -i " + inputVel + " -o " + outputDisp
    subprocess.call(svfExp + svfExpOptions, shell=True)

    #create warped Image
    warpImageOptions = " " + inputImage + " " + outputDisp + " " + outputImage
    subprocess.call(warpImage + warpImageOptions, shell=True)

    #Delete displacement field
    subprocess.call("rm " + outputDisp, shell=True)

    #Go to next file
    num += 1
    scale += options.stepSize

#combine images:
combineImagesOptions = " " + str(options.numOfScales) + " " + ".mha" + " " + resultsPath + "warpedByDisp" + " " + resultsPath + options.outputFile
subprocess.call(combineImages + combineImagesOptions, shell = True)
#delete inidvidual images:
subprocess.call("rm " + resultsPath + "warpedByDisp*.mha", shell=True)
