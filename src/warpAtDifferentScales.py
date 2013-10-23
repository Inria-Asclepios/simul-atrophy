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

#inputs of the executables
inputImage = dataPath + "Template/baselineDownsampleBy4.mha"
inputVel = resultsPath + "DownsampleBy4vel.mha"

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-n', '--number', dest='numOfScales', help='Total Number of images to produce', type=int)
parser.add_option('-i', '--initial', dest='startScale', help = 'Starting scale', type=int)
parser.add_option('-s', '--step', dest='stepSize', help = 'Step size between two scales', type=int)

(options, args) = parser.parse_args()

if options.numOfScales is None:
    options.numOfScales = int(raw_input('Enter total number of images to produce'))

if options.startScale is None:
    options.startScale = int(raw_input('Enter starting scale value:'))

if options.stepSize is None:
    options.stepSize = int(raw_input('Enter step size:'))


scale = options.startScale

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
combineImagesOptions = " " + str(options.numOfScales) + " " + ".mha" + " " + resultsPath + "warpedByDisp" + " " + resultsPath + "image4D.mha"
subprocess.call(combineImages + combineImagesOptions, shell = True)
#delete inidvidual images:
subprocess.call("rm " + resultsPath + "warpedByDisp*.mha", shell=True)
