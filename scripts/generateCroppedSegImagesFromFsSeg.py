#!/usr/bin/env python
import subprocess
import optparse
import os
import os.path as op


def print_and_execute(cmd):
    """
    Prints the command and executes it in the shell
    """
    print cmd + '\n'
    subprocess.call(cmd, shell=True)

#parse the command line inputs
parser = optparse.OptionParser()
parser.add_option('-m', '--mrOriginal', dest='mrOriginal', help='MRI T1 image. A cropped version of this image will be created ')
parser.add_option('-i', '--inImage', dest='inImage', help='input freesurfer aparc+aseg file in the same space as the input MR image.')
parser.add_option('-d', '--distThreshold', dest='distThreshold', help='distance from the tissue up to which  to fill CSF')
parser.add_option('-c', '--cropPadRadius', dest='cropPadRadius', help='Pad radius to use when cropping using the bounding box of the region that included newly added CSF regions.')
parser.add_option('-t', '--atrophyTables', dest='atrophyTables', help='a file that contains table with label values and the corresponding desired atrophy values')
parser.add_option('-r', '--resPath', dest='resPath', help='path to save results. Use . for current working directory.')
(options, args) = parser.parse_args()

AdLemModelPath = os.getenv('ADLEM_DIR')
antsPath = os.getenv('ANTS_BIN')
if AdLemModelPath is None:
    raise ValueError('environment variable ADLEM_DIR not found')
if antsPath is None:
    raise ValueError('environment variable ANTS_BIN not found.')

# Set required paths and binary aliases
# if options.isCluster is True:
#     antsPath = '/epi/asclepios2/bkhanal/antsbin/bin'
#     AdLemModelPath = '/epi/asclepios2/bkhanal/works/AdLemModel'
# else:
#     antsPath = '/home/bkhanal/Documents/softwares/antsbin/bin'
#     AdLemModelPath = '/home/bkhanal/works/AdLemModel'

if options.inImage is None:
    options.inImage = raw_input('Enter a valid aparc+aseg freesurfer segmentation')
inImage = options.inImage

if options.distThreshold is None:
    options.distThreshold = raw_input('distance from the tissue up to which  to fill CSF')
distThreshold = options.distThreshold

if options.resPath is None:
    options.resPath = raw_input('Enter directory where results will be saved; Use . to use current working directory')
resPath = options.resPath
if not op.exists(resPath):
    os.makedirs(resPath)

if options.cropPadRadius is None:
    options.cropPadRadius = raw_input('Enter pad radius for cropping. See help: -h for more info.')
cropPadRadius = options.cropPadRadius

if options.atrophyTables is None:
    options.atrophyTables = raw_input('Enter existing filename with a valid labels and corresponding atrophy values.')
atrophyTables = options.atrophyTables

if options.mrOriginal is None:
    options.mrOriginal = raw_input('Enter original T1 input structural MRI.')
t1Image = options.mrOriginal

# Output filenames
outMaskForModel = op.join(resPath, 'maskwithCsfD' + distThreshold + 'R' + cropPadRadius + '.nii.gz')
outLabelForAtrophy = op.join(resPath, 'labelsWithCsfD' + distThreshold + 'R' + cropPadRadius + '.nii.gz')
outAtrophyMap = op.join(resPath, 'atrophyMapD' + distThreshold + 'R' + cropPadRadius + '.nii.gz')
t1Cropped = op.join(resPath, 't1D' + distThreshold + 'R' + cropPadRadius + '.nii.gz')

# Binaries aliases
ImageMath = antsPath + '/ImageMath 3 '
binarize = AdLemModelPath + '/build/src/BinarizeThreshold -i '
getDistance = AdLemModelPath + '/build/src/signedDanielssonDistance '
crop = antsPath + '/ExtractRegionFromImageByMask 3 '

# Binarize input with 1. only tissue excluding all CSF to create distance map. 2. tissue with csf from freesurfer.
# With freesurfer csf
inBinaryWithCsf = op.join(resPath, 'tmp++InBinaryWithCsf.nii.gz')
command = binarize + inImage + ' -o ' + inBinaryWithCsf + ' -l 1 -u 15000'
print_and_execute(command)
# Now binary input excluding CSF (i.e. just tissue no CSF at all)
inBinary = op.join(resPath, 'tmp++InBinary.nii.gz')
# Extract freesurferCsf labels
freesurferCsf =op.join(resPath, 'tmp++freesurferCsf.nii.gz')
command = AdLemModelPath + '/build/src/createImageFromLabelImage -t ' + AdLemModelPath + '/configFiles/freeSurferCsfLabels -l ' \
          + inImage + ' -o ' + freesurferCsf
print_and_execute(command)
# Change freesurfer csf labels to prepare for binarization of input excluding freesurferCSF.
freesurferCsfMult = op.join(resPath, 'tmp++freesurferCsfMult.nii.gz')
command = ImageMath + freesurferCsfMult + ' m ' + freesurferCsf + ' -300 '  # 221 is max label for CSF region.
print_and_execute(command)
command = ImageMath + inBinary + ' + ' + inImage + ' ' + freesurferCsfMult
print_and_execute(command)
# Binarize by excluding freesurfer csf
command = binarize + inBinary + ' -o ' + inBinary + ' -l 1 -u 15000'
print_and_execute(command)

# get the distance of the binary input
distImage = op.join(resPath, 'tmp++distImage.nii.gz')
command = getDistance + inBinary + ' ' + distImage
print_and_execute(command)
# Extract desired CSF region input distance threshold
sulcalCsf = op.join(resPath, 'tmp++sulcalCsf.nii.gz')
command = binarize + distImage + ' -o ' + sulcalCsf + ' -l 0 -u ' + distThreshold
print_and_execute(command)

# Create a mask for the model i.e. labels as: NBR - 0  , CSF - 1, Tissue - 2
# Make tissue as 2
command = ImageMath + outMaskForModel + ' m ' + inBinary + ' 2'
print_and_execute(command)
# Include freesurfer csf regions
command = ImageMath + outMaskForModel + ' + ' + outMaskForModel + ' ' + freesurferCsf
print_and_execute(command)
# Now combine sulcal CSF regions with CSF label as 1 ( add only at those places where freesurfer had zero labels.)
command = ImageMath + outMaskForModel + ' addtozero ' + outMaskForModel + ' ' + sulcalCsf
print_and_execute(command)

# Combine CSF to the input image with CSF label AS 24.
command = ImageMath + sulcalCsf + ' m ' + sulcalCsf + ' 24 '
print_and_execute(command)
# Now combine sulcal CSF regions with CSF label as 1 ( add only at those places where freesurfer had zero labels.)
command = ImageMath + outLabelForAtrophy + ' addtozero ' + inImage + ' ' + sulcalCsf
print_and_execute(command)

# Now create atrophy map:
# The first atrophy table is used to create while the subsequent ones modify the
# one created with the first one.
commandPrefix = AdLemModelPath + '/build/src/createImageFromLabelImage -l ' + \
                outLabelForAtrophy + ' -o ' + outAtrophyMap + ' -t '
tables = atrophyTables.split()
# print('list of tables: ')
# print tables
command = commandPrefix + tables[0]
print_and_execute(command)
for table in tables[1:]: # from the remaining tables modify the output file
    command = commandPrefix + table + ' -m ' + outAtrophyMap
    print_and_execute(command)

# Crop all the images:
# Make a mask of including CSF regions with distance one more than the input distance.
cropMask = op.join(resPath, 'cropMaskD' + distThreshold + 'PadRadius' + cropPadRadius + '.nii.gz')
command = getDistance + inBinaryWithCsf + ' ' + distImage
print_and_execute(command)
command = binarize + distImage + ' -o ' + cropMask + ' -l -1000 -u ' + distThreshold
print_and_execute(command)

# Use this mask to crop the images.
command = crop + outAtrophyMap + ' ' + outAtrophyMap + ' ' + cropMask + ' 1 ' + cropPadRadius
print_and_execute(command)

command = crop + outMaskForModel + ' ' + outMaskForModel + ' ' + cropMask + ' 1 ' + cropPadRadius
print_and_execute(command)

command = crop + outLabelForAtrophy + ' ' + outLabelForAtrophy + ' ' + cropMask + ' 1 ' + cropPadRadius
print_and_execute(command)

command = crop + t1Image + ' ' + t1Cropped + ' ' + cropMask + ' 1 ' + cropPadRadius
print_and_execute(command)

# Delete all tmp++ prefixed files from the result directory.
command = 'rm ' + op.join(resPath, 'tmp++*')
print_and_execute(command)
