#!/usr/bin/env python
import subprocess
import sys
import optparse

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-i', '--inImage', dest='inImage', help='Filename of the input Segmentation whose convex hull is to be taken after optionally padding with CSF i.e. label=1')
parser.add_option('-m', '--inMask', dest='inMask', help='Filename of the input skull strip mask that was used in generating segementation mask')
parser.add_option('-o', '--outImage', dest='outImage', help='Filename of the output file')
parser.add_option('-w', '--width', dest='width', help='Voxel Width with which to dilate the input mask')

(options, args) = parser.parse_args()

if options.inImage is None:
    options.inImage = raw_input('Enter the input image FileName')

if options.inMask is None:
    options.inMask = raw_input('Enter the input Mask, e.g. a skull strip mask')

if options.outImage is None:
    options.outImage = raw_input('Enter the output image FileName')

if options.width is None:
    options.width = raw_input('Enter the mask dilation radius in voxel size')

#executables:
ImageMath = "ImageMath 3 "

outFileTmp = "dilated.nii.gz"
subprocess.call(ImageMath + outFileTmp + " MD " + options.inMask + " " + options.width, shell=True)
subprocess.call(ImageMath + outFileTmp + " - " + outFileTmp + " " + options.inMask, shell=True)
subprocess.call(ImageMath + outFileTmp + " + " + outFileTmp + " " + options.inImage, shell=True)
subprocess.call("ExtractRegionFromImage 3 " + outFileTmp + " " + options.outImage + " 1", shell=True)

print "now use: ExtractRegionFromImage 3 " + "inputImage outputImage indexDisplayedabove:axbxc size+index-1"
# $ ExtractRegionFromImage
# Usage 1: ExtractRegionFromImage ImageDimension inputImage outputImage minIndex maxIndex 
# Usage 2: ExtractRegionFromImage ImageDimension inputImage outputImage label 

# $ ExtractRegionFromImage 3 ~/data/patients/002_S_0938Atrophy1/002_S_0938_T0_W1_MNI_Stripped_segCsfPadded.nii.gz ~/data/patients/002_S_0938Atrophy1/002_S_0938_T0_W1_MNI_Stripped_segCsfPaddedConvexHull.nii.gz 1
# ImageRegion (0x7fffcfce0300)
#   Dimension: 3
#   Index: [17, 16, 6]
#   Size: [147, 187, 151]

# $ ExtractRegionFromImage 3 ~/data/patients/002_S_0938Atrophy1/002_S_0938_T0_W1_MNI_Stripped.nii.gz ~/data/patients/002_S_0938Atrophy1/002_S_0938_T0_W1_MNI_Stripped_ConvexHull.nii.gz 17x16x6 163x202x156
# ImageRegion (0x7fffae1c6a60)
#   Dimension: 3
#   Index: [17, 16, 6]
#   Size: [147, 187, 151]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Subsampled case:
# $ExtractRegionFromImage 3 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb_ConvexHull.nii.gz 8x8x3 82x101x78
# ImageRegion (0x7fff96816760)
#   Dimension: 3
#   Index: [8, 8, 3]
#   Size: [75, 94, 76]

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Creating individual masks of GM, WM and computing volume fractions and assigning atrophy maps accordingly:
# $ fslmaths -dt int tstMask.nii.gz -thr 3 tstMaskWM.nii.gz
# $ fslmaths -dt int tstMaskWM.nii.gz -div 3 tstMaskWM.nii.gz

# $ fslmaths -dt int tstMask.nii.gz -uthr 2 tstMaskGmTmp.nii.gz
# $ fslmaths -dt int tstMaskGmTmp.nii.gz -thr 2 tstMaskGM.nii.gz
# $ fslmaths -dt int tstMaskGM.nii.gz -div 2 tstMaskGM.nii.gz

# $ fslstats tstMaskGM.nii.gz -V
# 85150 681200.000000 
# $ fslstats tstMaskWM.nii.gz -V
# 72868 582944.000000 


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Adapting the cerebellum mask:
# $ fslmaths ../results/patients/002_S_0938/002_S_0938_Cereb.nii.gz -subsamp2 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -thr 0.5 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -sub 1 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -abs ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
# $ fslmaths ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz -thr 0.9 ../results/patients/002_S_0938SubSampled2/002_S_0938SubSampled2_Cereb.nii.gz
