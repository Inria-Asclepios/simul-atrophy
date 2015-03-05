#!/usr/bin/env python
import subprocess
import sys
import optparse

def printAndExecute(command):  #prints the command and executes it in the shell
    print command + '\n'
    subprocess.call(command, shell=True)

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('--isCluster', dest='isCluster', action = 'store_true', help='given: use nef cluster paths. not given: use laptop paths')
parser.add_option('-p', '--filePrefix', dest='filePrefix', help='prefix to be added to all the output files.')
parser.add_option('-i', '--inImage', dest='inImage', help='input T1 structural MRI.')
parser.add_option('--noSkullStripping', dest='noSkullStripping', action = 'store_true', help='given: skip skull stripping and use the skull strip mask given using the option -s. Not given: do skull stripping')
parser.add_option('-s', '--ssMask', dest='ssMask', help='Binary skull strip mask defining the brain region of the input structural MRI with values 1.')
parser.add_option('--runFast', dest='runFast', action = 'store_true', help='given: run fsl fast to segment gm/wm/csf. Not given: uses the file prefix given to infer the fast segmented files')
parser.add_option('-H', '--segSmoothness', dest='segSmoothness', help='fsl fast option -H that determines smoothness of the segmented GM/WM/CSF. Default=0.3; Must be provided to determine fast result prefix.')
parser.add_option('--runFirstAll', dest='runFirstAll', action = 'store_true', help='given: run fsl run_first_all to segment sub-cortical structures. Not given: uses the file prefix given to infer the first segmented files')
parser.add_option('-w', '--dilateWidth', dest='dilateWidth', help='Voxel Width with which to dilate the brain mask which will determine the width with which to surround the segmented image by CSF; Set 0 if you dont want to dilate')
parser.add_option('-c', '--beforeCropPadWidth', dest='beforeCropPadWidth', help=' Width which will be added to the extracted bounding box region before cropping the input image.')

(options, args) = parser.parse_args()

if options.filePrefix is None:
    options.filePrefix = raw_input('Enter a prefix to be added to all the result files.')
filePrefix = options.filePrefix

if options.inImage is None:
    options.inImage = raw_input('Enter a valid structural MRI.')
inImage = options.inImage

if options.dilateWidth is None:
    options.dilateWidth = raw_input('Enter the mask dilation radius in voxel size. 0 if you dont want to dilate')

if options.beforeCropPadWidth is None:
    options.beforeCropPadWidth = raw_input('Enter the pad radius to be used before cropping.')

# Set the paths:
if options.isCluster is True:
    skullStrip = '/epi/asclepios2/bkhanal/ROBEX/runROBEX.sh '
    antsPath = '/epi/asclepios2/bkhanal/antsbin/bin'
    fslPath = '/epi/asclepios2/bkhanal/fsl/bin'
    AdLemModelPath = '/epi/asclepios2/bkhanal/works/AdLemModel'
    myToolsPath = '/epi/asclepios2/bkhanal/works/tools'
else:
    skullStrip = '/home/bkhanal/Documents/softwares/ROBEX/runROBEX.sh '
    antsPath = '/home/bkhanal/Documents/softwares/antsbin/bin'
    fslPath = '/home/bkhanal/Documents/softwares/fsl/fsl/bin'
    AdLemModelPath = '/home/bkhanal/works/AdLemModel'
    myToolsPath = '/home/bkhanal/works/tools'


# Command aliases
binarize = AdLemModelPath + '/build/src/BinarizeThreshold -i '
ImageMath = antsPath + '/ImageMath 3 '

# If asked do skull stripping, otherwise ask for a brain mask.
ssImage = filePrefix + '_ss.nii.gz'
if options.noSkullStripping is True:
    if options.ssMask is None:
        options.ssMask = raw_input('Enter a valid skull strip mask that defines brain region of the input structural MRI')
    ssMask = options.ssMask
    command = fslPath + '/fslmaths ' + inImage + ' -mas ' + ssMask + ' ' + ssImage
    printAndExecute(command)
else:
    # Skull strip mask
    ssMask = filePrefix + '_ssMask.nii.gz'
    command = skullStrip + inImage + ' ' + ssImage + ' ' + ssMask
    printAndExecute(command)

# If asked, do segmentations using FAST on the skull stripped images:
fastOutputPrefix = filePrefix
if options.runFast is True:
    if options.segSmoothness is None:
        options.segSmoothness = '0.3'
    command = fslPath + '/fast -g --nopve -H ' + options.segSmoothness + ' -o ' + fastOutputPrefix + ' ' + ssImage
    printAndExecute(command)

# If asked, do cortical substructures segmentation of ssImage
if options.runFirstAll is True:
    command = fslPath + '/run_first_all -m fast -b -i ' + ssImage + ' -o ' + filePrefix
    printAndExecute(command)
csSeg = filePrefix + '_all_fast_firstseg.nii.gz'


#If label table is given, use it to create atrophy map.
# Replace label code and atrophy generation code.


# Create the segMask for model i.e. with labels: NBR -> 0, CSF -> 1, GM/WM -> 2:
fastSeg = fastOutputPrefix + '_seg.nii.gz'
modelSegTmp = 'tmp++segMaskForModel.nii.gz'
printAndExecute(binarize + fastSeg + ' -o ' + modelSegTmp + ' -l 2 -u 3 -x 2 -s true')

# Fill CSF around the segmentation if non zero mask dilation width is given
if (options.dilateWidth == '0'):
    cropMask = ssMask
else:
    # Dilate ssMask which will be used to crop segmented images.
    ssMaskDilated = filePrefix + '_ssMaskDilated.nii.gz'
    # Dilate the input mask with the given radius.
    printAndExecute(ImageMath + ssMaskDilated + " MD " + ssMask + " " + options.dilateWidth)
    #Substract the input mask from the dilated image to get a mask of newly added pixel locations.
    outFileTmp = 'tmp++Dilation.nii.gz'
    printAndExecute(ImageMath + outFileTmp + " - " + ssMaskDilated + " " + ssMask)
    #Add the result to the segmented mask image for the model. This means extra CSF gets added at the newly created locations due to dilation of the brain mask.
    printAndExecute(ImageMath + modelSegTmp + " + " + outFileTmp + " " + modelSegTmp)
    cropMask = ssMaskDilated


#Extract convex hull bounding box (BB) of all the images that are required by the model.
crop = antsPath + '/ExtractRegionFromImageByMask 3 '
cropPostfix = 'CW' + options.dilateWidth + '.nii.gz'

inImageCropped = filePrefix + cropPostfix
printAndExecute(crop + inImage + " " + inImageCropped + " " + cropMask + " 1 " + options.beforeCropPadWidth)

modelSeg = filePrefix + '_segMaskForModelCW' + options.dilateWidth + '.nii.gz'
printAndExecute(crop + modelSegTmp + " " + modelSeg + " " + cropMask + " 1 " + options.beforeCropPadWidth)

# Delete all tmp++ prefixed files from the result directory.
printAndExecute('rm ' + 'tmp++*')

