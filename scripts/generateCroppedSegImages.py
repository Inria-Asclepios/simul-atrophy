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


(options, args) = parser.parse_args()

if options.filePrefix is None:
    options.filePrefix = raw_input('Enter a prefix to be added to all the result files.')
filePrefix = options.filePrefix

if options.inImage is None:
    options.inImage = raw_input('Enter a valid structural MRI.')
inImage = options.inImage


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

ssImage = filePrefix + '_ss.nii.gz'
if options.noSkullStripping is True:
    if options.ssMask is None:
        options.ssMask = raw_input('Enter a valid skull strip mask that defines brain region of the input structural MRI')
    ssMask = options.ssMask
    command = fslPath + '/fslmaths ' + inImage + ' -mas ' + ssMask + ' ' + ssImage
    #fslmaths ../results/patients/testGenerateScript/T1.nii.gz -mas ../results/patients/testGenerateScript/T1_ssMask.nii.gz ../results/patients/testGenerateScript/testingSS.nii.gz
    printAndExecute(command)
else:
    # Skull strip mask
    ssMask = filePrefix + '_ssMask.nii.gz'
    command = skullStrip + inImage + ' ' + ssImage + ' ' + ssMask
    #  ./ROBEX ~/works/AdLemModel/results/patients/003_S_4288/T1.nii.gz ~/works/AdLemModel/results/patients/003_S_4288/T1_ss.nii.gz ~/works/AdLemModel/results/patients/003_S_4288/T1_ssMask.nii.gz
    printAndExecute(command)


# Do the required segmentations on the skull stripped images:
# GM, WM, CSF segmentation:
if options.runFast is True:
    if options.segSmoothness is None:
        options.segSmoothness = '0.3'
fastOutputPrefix = filePrefix + 'H' + options.segSmoothness
command = fslPath + '/fast -g --nopve -H ' + options.segSmoothness + ' -o ' + fastOutputPrefix + ' ' + ssImage
#   863  fast -g --nopve -H 0.5 -o T1H0_5 T1_ss.nii.gz
printAndExecute(command)
csfSeg = fastOutputPrefix + '_seg_0.nii.gz'
gmSeg = fastOutputPrefix + '_seg_1.nii.gz'
wmSeg = fastOutputPrefix + '_seg_2.nii.gz'

# Cortical substructures segmentation of ssImage except brainStem
command = fslPath + '/run_first_all -m fast -b -i ' + ssImage + ' -o ' + filePrefix
#  1011  run_first_all -m fast -b -s L_Accu,L_Amyg,L_Caud,L_Hipp,L_Pall,L_Puta,L_Thal,R_Accu,R_Amyg,R_Caud,R_Hipp,R_Pall,R_Puta,R_Thal -i ../results/patients/003_S_4288/inputProcessing/T1_ss.nii.gz -o ../results/patients/003_S_4288/inputProcessing/T1H0_5noStemAll
printAndExecute(command)
csSeg = filePrefix + '_all_fast_firstseg.nii.gz'

#Separate out the stem and nonStem parts and binarize them:
# Non-stem part
csNoStemSeg = filePrefix + '_all_fast_firstNoStemSeg.nii.gz'
command = AdLemModelPath + '/build/src/BinarizeThreshold -i ' + csSeg + ' -o ' + csNoStemSeg + ' -l 16 -u 16.5 -x 0 -s true'
printAndExecute(command)
# 1102  src/BinarizeThreshold -i ../results/patients/003_S_4288TestFsl/T1_all_fast_firstseg.nii.gz -o ../results/patients/003_S_4288TestFsl/T1_all_fast_firstsegNoStem.nii.gz -l 16 -u 16.5 -x 0 -s true
command = AdLemModelPath + '/build/src/BinarizeThreshold -i ' + csNoStemSeg + ' -o ' + csNoStemSeg + ' -l 1 -u 100 -x 1'
# 1103  src/BinarizeThreshold -i ../results/patients/003_S_4288TestFsl/T1_all_fast_firstsegNoStem.nii.gz -o ../results/patients/003_S_4288TestFsl/T1_all_fast_firstsegNoStemBinary.nii.gz -l 1 -u 100 -x 1
printAndExecute(command)

#Stem part
stemSeg = filePrefix + '_all_fast_firstStemSeg.nii.gz'
command = AdLemModelPath + '/build/src/BinarizeThreshold -i ' + csSeg + ' -o ' + stemSeg + ' -l 16 -u 16.5 -x 1'
# src/BinarizeThreshold -i ../results/patients/003_S_4288TestFsl/T1_all_fast_firstseg.nii.gz -o ../results/patients/003_S_4288TestFsl/T1_all_fast_firstsegStemBinary.nii.gz -l 16 -u 16.5 -x 1
printAndExecute(command)


# Make a binary mask with object as GM+WM+CS and everything else including CSF as the background.
gmWmCsNoStemSeg = fastOutputPrefix + '_seg_12_firstAllNoStem.nii.gz' # 3 for CS without stem. Use 4 for CS with stem.

# First combine the cs structures except stem and gm and wm which is still a binary image. 
command = antsPath + '/ImageMath 3 ' + gmWmCsNoStemSeg + ' + ' +  gmSeg + ' ' + wmSeg
printAndExecute(command)
command = antsPath + '/ImageMath 3 ' + gmWmCsNoStemSeg + ' addtozero ' +  gmWmCsNoStemSeg + ' ' + csNoStemSeg
printAndExecute(command)

#Now remove the brain stem
tmpImage = 'tmp.nii.gz'
command = antsPath + '/ImageMath 3 ' + tmpImage + ' m ' +  stemSeg + ' 5 '
printAndExecute(command)
command = antsPath + '/ImageMath 3 ' + gmWmCsNoStemSeg + ' + ' +  gmWmCsNoStemSeg + ' ' + tmpImage
printAndExecute(command)
command = AdLemModelPath + '/build/src/BinarizeThreshold -i ' + gmWmCsNoStemSeg + ' -o ' + gmWmCsNoStemSeg + ' -l 5 -u 7 -x 0 -s true'
printAndExecute(command)

command = 'rm ' + tmpImage
printAndExecute(command)

# If separate atrophy values are desired at different segmentations, use the individual segmentations that are saved and multiply them before adding such as:
# Note: Do all these before cropping because that way, the cropping and padding can be done at the very end on only those images that are required by the AdLemModel code!
# ;; Atrophy generation:
# ;; Combine gmSeg with corticalSturcturesSeg excluding the brain stem and add desired atrophy.
#  1019  ImageMath 3 tmp.nii.gz m ../results/patients/003_S_4288/inputProcessing/T1H0_5_segP4C0_1.nii.gz -0.05 ;;GM atrophy 
#  1020  ImageMath 3 tmp1.nii.gz m ../results/patients/003_S_4288/inputProcessing/T1H0_5noStemAll_all_fast_firstsegBinaryP4C0.nii.gz -0.1  ;;cs atrophy
#  1022  ImageMath 3 ../results/patients/003_S_4288GmAtrophy/atrophyGm5pCstructs10p.nii.gz + tmp.nii.gz tmp1.nii.gz  ;;combine the atrophies.


# ;; White matter distance from boundary in negative. All other zero; Set the brain stem to zero too.
distImage = fastOutputPrefix + 'DistMap.nii.gz'
command = '/home/bkhanal/works/quickTests/build/src/quickTest ' + wmSeg + ' ' + distImage
#  897  src/quickTest ~/works/AdLemModel/results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2.nii.gz ~/works/AdLemModel/results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImage.nii.gz 
printAndExecute(command)
command = antsPath + '/ImageMath 3 ' + distImage + ' - ' + distImage + ' 1' # make all WM region negative
#  1047  ImageMath 3 ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImageClippedMinus1.nii.gz - ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImageClipped.nii.gz 1
printAndExecute(command)
command = AdLemModelPath + '/build/src/BinarizeThreshold -i ' + distImage + ' -o ' + distImage + ' -l 0 -u 1000 -x 0 -s true'  # set all non-WM region to zero.
#  1048  src/BinarizeThreshold -i ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImageClippedMinus1.nii.gz -o ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImageClippedThresholded.nii.gz -l 0 -u 10 -x 0 -s true 
printAndExecute(command)

# Set brain stem area to zero.
tmpImage = 'tmpImage.nii.gz'
command = antsPath + '/ImageMath 3 ' + tmpImage + ' m ' + stemSeg + ' 1000'
# 1058  ImageMath 3 tmp.nii.gz m ../results/patients/003_S_4288/inputProcessing/T1H0_5brainStemP4C0.nii.gz 10
printAndExecute(command)
command = antsPath + '/ImageMath 3 ' + distImage + ' + ' + tmpImage + ' ' + distImage
# 1062  ImageMath 3 tmp1.nii.gz addtozero tmp.nii.gz ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImageClippedThresholdedP4C0.nii.gz 
printAndExecute(command)
command = AdLemModelPath + '/build/src/BinarizeThreshold -i ' + distImage + ' -o ' + distImage + ' -l 0 -u 1000 -x 0 -s true'
# src/BinarizeThreshold -i tmp1.nii.gz -o ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceMapClippedThresholdedNobrainStemP4C0.nii.gz -l 1 -u 200 -x 0 -s true
printAndExecute(command)

command = 'rm ' + tmpImage
printAndExecute(command)

# Dilate ssMask and which will be used to crop segmented images.
ssMaskDilated = filePrefix + '_ssMaskDilated.nii.gz'

# ;;pad and crop the ssImage. First add CSF around segmentedImage with the given csfWidth (=4 here). This is done with dilation of the ssMask. Then extract bounding box using the ssMaskDilated after again padding with given cropPad width(=0 here).
#   872  ../scripts/padCsfInStrippedSkullAndTakeConvexHull.py -i ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg.nii.gz -m ../results/patients/003_S_4288/inputProcessing/T1_ssMask.nii.gz -w 4 -v 1 -p 0 -o ../results/patients/003_S_4288/inputProcessing/T1H0_5_segP4C0.nii.gz -d ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz

# ;; Crop all the required images (for atrophy generation and as model inputs) using ssMaskDilated.
#  1007  ExtractRegionFromImageByMask 3 ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_0.nii.gz ../results/patients/003_S_4288/inputProcessing/T1H0_5_segP4C0_0.nii.gz ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz 1 0
#  1008  ExtractRegionFromImageByMask 3 ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_1.nii.gz ../results/patients/003_S_4288/inputProcessing/T1H0_5_segP4C0_1.nii.gz ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz 1 0
#  1009  ExtractRegionFromImageByMask 3 ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2.nii.gz ../results/patients/003_S_4288/inputProcessing/T1H0_5_segP4C0_2.nii.gz ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz 1 0
#  1016  ExtractRegionFromImageByMask 3 ../results/patients/003_S_4288/inputProcessing/T1H0_5noStemAll_all_fast_firstseg.nii.gz ../results/patients/003_S_4288/inputProcessing/T1H0_5noStemAll_all_fast_firstsegP4C0.nii.gz ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz 1 0
#  1027  ExtractRegionFromImageByMask 3 ../results/patients/003_S_4288/inputProcessing/T1.nii.gz ../results/patients/003_S_4288/inputProcessing/T1P4C0.nii.gz ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz 1 0
#  1061  ExtractRegionFromImageByMask 3 ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImageClippedThresholded.nii.gz ../results/patients/003_S_4288/inputProcessing/T1H0_5_seg_2DistanceImageClippedThresholdedP4C0.nii.gz ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz 1 0                                                                                                                                                                                                                                                     
#  1030  src/ExtractRegionOfDtiByMask ../results/patients/003_S_4288/inputProcessing/dt.nii.gz ../results/patients/003_S_4288/inputProcessing/dtP4C0.nii.gz ../results/patients/003_S_4288/inputProcessing/T1_ssMaskDilatedW4.nii.gz 1 0
#  1032  src/SwitchUpperLowerTriangular ../results/patients/003_S_4288/inputProcessing/dtP4C0.nii.gz ../results/patients/003_S_4288/inputProcessing/dt_itkP4C0.nii.gz


