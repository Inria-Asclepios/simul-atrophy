#!/usr/bin/env python
import subprocess
import sys
import optparse

#Executables that must exist:
#ImageMath
sliceFrom3D = "/home/bkhanal/works/tools/marcoSliceExtractor/myImgSliceExtractor"

#convert (from ImageMagick)
#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-i', '--inputFile', dest='inputFile', help='Filename of the input 4D image')
parser.add_option('-o', '--outputFile', dest='outputFile', help='Filename of the output file')
parser.add_option('-n', '--numOfTimeSteps', dest='numOfTimeSteps', help='Number of time steps present in the input file',type='int')
parser.add_option('-a', '--sliceAxis', dest='sliceAxis', help='Slice axis: 0,1 or 2')
parser.add_option('-s', '--sliceNumber', dest='sliceNumber', help='slice number to be extracted')
parser.add_option('-r', '--resize', dest='resize', help='resize to maximum possible size in each direction, e.g. 400x400')
parser.add_option('-d', '--delay', dest='delay', help='time delay between frames in milliseconds')

(options, args) = parser.parse_args()

if options.inputFile is None:
    options.inputFile = raw_input('Enter the input image FileName')

if options.outputFile is None:
    options.outputFile = raw_input('Enter the output image FileName')

if options.numOfTimeSteps is None:
    options.numOfTimeSteps = int(raw_input('Enter number of times steps'))

if options.sliceAxis is None:
    options.sliceAxis = raw_input('Enter slice axis 0,1 or 2')

if options.sliceNumber is None:
    options.sliceNumber = raw_input('Enter slice number to be extracted')

if options.resize is None:
    options.resize = raw_input('Enter size of output gif file; e.g. 400x400')

if options.delay is None:
    options.delay = raw_input('Enter dealy in frames in output gif file; in ms e.g. 10')


#split 4d image to get all the 3D images:
subprocess.call("ImageMath 4 tmp3DImage.nii.gz TimeSeriesDisassemble " + options.inputFile, shell=True)

axis = options.sliceAxis
num = 0
print 'number of time steps: ' + str(options.numOfTimeSteps) + '\n'
while num < options.numOfTimeSteps:
    #outputs/results of the executables
    index = str(num+100) #This is because ImageMath saves individual slices starting from 100.
    tmp3DImage = 'tmp3DImage' + index + '.nii.gz'
    tmp2DImage = 'slice' + index + '.png'
    extractSlice = sliceFrom3D + " " + tmp3DImage + " " + axis + " " + str(options.sliceNumber) + " " + tmp2DImage
    #print 'command: ' + extractSlice + '\n'
    subprocess.call(extractSlice, shell=True)
    subprocess.call("convert -rotate 180 " + tmp2DImage + " " + tmp2DImage, shell=True)

    #Delete individual 3D files.
    deleteTmpFiles = "rm tmp3DImage" + index + ".nii.gz"
    subprocess.call(deleteTmpFiles, shell=True)
    #Go to next file
    num += 1

#Now make the animation and delete individual 2D slices:
subprocess.call('convert slice1*.png -resize ' + options.resize + ' -set delay ' + options.delay + ' ' + options.outputFile, shell=True)
subprocess.call('rm slice1*.png', shell=True)


    


