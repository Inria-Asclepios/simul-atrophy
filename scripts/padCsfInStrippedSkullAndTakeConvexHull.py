#!/usr/bin/env python
import subprocess
import sys
import optparse

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-i', '--inImage', dest='inImage', help='Filename of the input Segmentation whose convex hull is to be taken after optionally padding with desired pixel value')
parser.add_option('-m', '--inMask', dest='inMask', help='Filename of the input skull strip mask that was used in generating segementation mask')
parser.add_option('-w', '--dilateWidth', dest='dilateWidth', help='Voxel Width with which to dilate the input mask which will determine the width with which to surround the input image by adding input pixelValue')
parser.add_option('-v', '--pixelAddValue', dest='pixelAddValue', help='Default = 1. Pixel value that will be added at all newly created non-zero location due to dilation of the input mask.')
parser.add_option('-p', '--beforeCropPadWidth', dest='beforeCropPadWidth', help=' Defalut = 0. Width which will be added to the extracted bounding box region before cropping the input image.')
parser.add_option('-o', '--outImage', dest='outImage', help='Filename of the output file')
parser.add_option('-d', '--dilatedMask', dest='dilatedMask', help='Filename of the output dilated mask, that could be used to extract regions for other images to have same size as the one extracted with this script.')

(options, args) = parser.parse_args()

if options.inImage is None:
    options.inImage = raw_input('Enter the input image FileName')

if options.inMask is None:
    options.inMask = raw_input('Enter the input Mask, e.g. a skull strip mask')

if options.dilateWidth is None:
    options.dilateWidth = raw_input('Enter the mask dilation radius in voxel size')

if options.pixelAddValue is None:
    options.pixelAddValue = 1

if options.beforeCropPadWidth is None:
    options.beforeCropPadWidth = 0

if options.outImage is None:
    options.outImage = raw_input('Enter the output image fileName')

if options.dilatedMask is None:
    options.dilatedMask = raw_input('Enter the output dilated image mask fileName')


#executables:
ImageMath = "ImageMath 3 "

outFileTmp = "tmp.nii.gz"
dilatedMask = options.dilatedMask
# Assumptions:
# 1. The input mask is a binary with 1 being forground pixel.

# Dilate the input mask with the given radius.
subprocess.call(ImageMath + dilatedMask + " MD " + options.inMask + " " + options.dilateWidth, shell=True)

#Substract the input mask from the dilated image to get a mask of newly added pixel locations.
subprocess.call(ImageMath + outFileTmp + " - " + dilatedMask + " " + options.inMask, shell=True)

#Multiply the result by the pixel Value to be added to the input image at these locations.
subprocess.call(ImageMath + outFileTmp + " m " + outFileTmp + " " + options.pixelAddValue, shell=True)

#Add the result to the input image. This means pixelValue gets added at the newly created locations due to dilation of the input mask.
#This for example can pad a segmented brain image with csf around the skull stripped boundary.
subprocess.call(ImageMath + outFileTmp + " + " + outFileTmp + " " + options.inImage, shell=True)

#Extract convex hull bounding box (BB) of the input image using dilated mask after adding the input crop width to the BB.
subprocess.call("ExtractRegionFromImageByMask 3 " + outFileTmp + " " + options.outImage + " " + dilatedMask + " 1 " + options.beforeCropPadWidth, shell=True)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Adapting the cerebellum mask:
# $ fslmaths ../results/patients/002_S_0938/002_S_0938_Cereb.nii.gz -subsamp2 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -thr 0.5 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -sub 1 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -abs ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -thr 0.9 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
