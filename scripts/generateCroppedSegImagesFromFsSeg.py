#!/usr/bin/env python
import subprocess
import optparse

def printAndExecute(command):  #prints the command and executes it in the shell
    print command + '\n'
    subprocess.call(command, shell=True)

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('--isCluster', dest='isCluster', action = 'store_true', help='given: use nef cluster paths. not given: use laptop paths')
parser.add_option('-m', '--mrOriginal', dest='mrOriginal', help='MRI T1 image. A cropped version of this image will be created ')
parser.add_option('-i', '--inImage', dest='inImage', help='input freesurfer aparc+aseg file in the same space as the input MR image.')
parser.add_option('-d', '--distThreshold', dest='distThreshold', help='distance from the tissue up to which  to fill CSF')
parser.add_option('-c', '--cropPadRadius', dest='cropPadRadius', help='Pad radius to use when cropping using the bounding box of the region that included newly added CSF regions.')
parser.add_option('-t', '--atrophyTable', dest='atrophyTable', help='a file that contains table with label values and the corresponding desired atrophy values')
parser.add_option('-r', '--resPath', dest='resPath', help='path to save results. Use . for current working directory.')
(options, args) = parser.parse_args()

# Set required paths and binary aliases
if options.isCluster is True:
    #runRobex = '/epi/asclepios2/bkhanal/ROBEX/runROBEX.sh '
    antsPath = '/epi/asclepios2/bkhanal/antsbin/bin'
    #fslPath = '/epi/asclepios2/bkhanal/fsl/bin'
    AdLemModelPath = '/epi/asclepios2/bkhanal/works/AdLemModel'
    #myToolsPath = '/epi/asclepios2/bkhanal/works/tools'
    signedMaurerDistance = '/epi/asclepios2/bkhanal/works/tools/...'
else:
    #runRobex = '/home/bkhanal/Documents/softwares/ROBEX/runROBEX.sh '
    antsPath = '/home/bkhanal/Documents/softwares/antsbin/bin'
    #fslPath = '/home/bkhanal/Documents/softwares/fsl/fsl/bin'
    AdLemModelPath = '/home/bkhanal/works/AdLemModel'
    #myToolsPath = '/home/bkhanal/works/tools'
    getDistance = '/home/bkhanal/works/quickTests/build/src/signedMaurerDistance '


if options.inImage is None:
    options.inImage = raw_input('Enter a valid aparc+aseg freesurfer segmentation')
inImage = options.inImage

if options.distThreshold is None:
    options.distThreshold = raw_input('distance from the tissue up to which  to fill CSF')
distThreshold = options.distThreshold

if options.resPath is None:
    options.resPath = raw_input('Enter directory where results will be saved; Use . to use current working directory')
resPath = options.resPath.rstrip('/')  # I'll add slash if present, so remove it if given by the user.

if options.cropPadRadius is None:
    options.cropPadRadius = raw_input('Enter pad radius for cropping. See help: -h for more info.')
cropPadRadius = options.cropPadRadius

if options.atrophyTable is None:
    options.atrophyTable = raw_input('Enter existing filename with a valid labels and corresponding atrophy values.')
atrophyTable = options.atrophyTable

if options.mrOriginal is None:
    options.mrOriginal = raw_input('Enter original T1 input structural MRI.')
t1Image = options.mrOriginal

# Output filenames
outMaskForModel = resPath + '/maskwithCsfD' + distThreshold + 'R' + cropPadRadius + '.nii.gz'
outLabelForAtrophy = resPath + '/labelsWithCsfD' + distThreshold + 'R' + cropPadRadius + '.nii.gz'
outAtrophyMap = resPath + '/atrophyMapD' + distThreshold + 'R' + cropPadRadius + '.nii.gz'
t1Cropped = resPath + '/t1D' + distThreshold + 'R' + cropPadRadius + '.nii.gz'

# Binaries aliases
ImageMath = antsPath + '/ImageMath 3 '
binarize = AdLemModelPath + '/build/src/BinarizeThreshold -i '


# Binarize input with 1. only tissue excluding all CSF to create distance map. 2. tissue with csf from freesurfer.
# With freesurfer csf
inBinaryWithCsf = 'tmp++InBinaryWithCsf.nii.gz'
command = binarize + inImage + ' -o ' + inBinaryWithCsf + ' -l 1 -u 15000'
printAndExecute(command)
# Now binary input excluding CSF (i.e. just tissue no CSF at all)
inBinary = 'tmp++InBinary.nii.gz'
# Extract freesurferCsf labels
freesurferCsf ='tmp++freesurferCsf.nii.gz'
command = AdLemModelPath + '/build/src/replaceLabelsInLabelImage ' + AdLemModelPath + '/configFiles/freeSurferCsfLabels ' \
          + inImage + ' ' + freesurferCsf
printAndExecute(command)
# Change freesurfer csf labels to prepare for binarization of input excluding freesurferCSF.
freesurferCsfMult = 'tmp++freesurferCsfMult.nii.gz'
command = ImageMath + freesurferCsfMult + ' m ' + freesurferCsf + ' -300 '  # 221 is max label for CSF region.
printAndExecute(command)
command = ImageMath + inBinary + ' + ' + inImage + ' ' + freesurferCsfMult
printAndExecute(command)
# Binarize by excluding freesurfer csf
command = binarize + inBinary + ' -o ' + inBinary + ' -l 1 -u 15000'
printAndExecute(command)

# get the distance of the binary input
distImage = 'tmp++distImage.nii.gz'
command = getDistance + inBinary + ' ' + distImage
printAndExecute(command)
# Extract desired CSF region input distance threshold
sulcalCsf = 'tmp++sulcalCsf.nii.gz'
command = binarize + distImage + ' -o ' + sulcalCsf + ' -l 1 -u ' + distThreshold
printAndExecute(command)

# Create a mask for the model i.e. labels as: NBR - 0  , CSF - 1, Tissue - 2
# Make tissue as 2
command = ImageMath + outMaskForModel + ' m ' + inBinary + ' 2'
printAndExecute(command)
# Include freesurfer csf regions
command = ImageMath + outMaskForModel + ' + ' + outMaskForModel + ' ' + freesurferCsf
printAndExecute(command)
# Now combine sulcal CSF regions with CSF label as 1 ( add only at those places where freesurfer had zero labels.)
command = ImageMath + outMaskForModel + ' addtozero ' + outMaskForModel + ' ' + sulcalCsf
printAndExecute(command)

# Combine CSF to the input image with CSF label AS 24.
command = ImageMath + sulcalCsf + ' m ' + sulcalCsf + ' 24 '
printAndExecute(command)
# Now combine sulcal CSF regions with CSF label as 1 ( add only at those places where freesurfer had zero labels.)
command = ImageMath + outLabelForAtrophy + ' addtozero ' + inImage + ' ' + sulcalCsf
printAndExecute(command)

# Now create atrophy map:
command = AdLemModelPath + '/build/src/replaceLabelsInLabelImage ' + atrophyTable + ' ' + outLabelForAtrophy + ' ' + \
          outAtrophyMap
printAndExecute(command)

# Crop all the images:
# Make a mask of including CSF regions with distance one more than the input distance.
cropMask = resPath + '/cropMaskD' + distThreshold + 'PadRadius' + cropPadRadius + '.nii.gz'
command = getDistance + inBinaryWithCsf + ' ' + distImage
printAndExecute(command)
command = binarize + distImage + ' -o ' + cropMask + ' -l -1000 -u ' + distThreshold
printAndExecute(command)

# Use this mask to crop the images.
crop = antsPath + '/ExtractRegionFromImageByMask 3 '

command = crop + outAtrophyMap + ' ' + outAtrophyMap + ' ' + cropMask + ' 1 ' + cropPadRadius
printAndExecute(command)

command = crop + outMaskForModel + ' ' + outMaskForModel + ' ' + cropMask + ' 1 ' + cropPadRadius
printAndExecute(command)

command = crop + outLabelForAtrophy + ' ' + outLabelForAtrophy + ' ' + cropMask + ' 1 ' + cropPadRadius
printAndExecute(command)

command = crop + t1Image + ' ' + t1Cropped + ' ' + cropMask + ' 1 ' + cropPadRadius
printAndExecute(command)

# Delete all tmp++ prefixed files from the result directory.
command = 'rm ' + 'tmp++*'
printAndExecute(command)
