#!/usr/bin/python
"""
generateCroppedSegImagesFromFsSeg
Generate segmented images and masks for Adlem model.
"""
# Author: Bishesh Khanal <bishesh.khanal@inria.fr>
#
#

import argparse as ag
import os
import os.path as op
import bish_utils as bu

# Paths of required binaries
adlem_dir, ants_dir = None, None
# Required binaries
img_math, binarize, get_dist = None, None, None
img_from_label_img = None

def get_input_options():
    """Parse command line inputs """
    parser = ag.ArgumentParser(
        'Generate segmented images, mask and atrophy map that are adapted to be'
        ' inputs for AdLemModel.')
    parser.add_argument('t1_img', help='MRI T1 image. A cropped version of'
                        ' this image will be created ')
    parser.add_argument('in_seg', help='input freesurfer segmentation file in '
                        'the same space as the input MR image.')
    parser.add_argument(
        'fs_csf_labels', help='File with FS labels to be considered as'
        ' CSF regions for the model. Format: 1st columns contain labels, second'
        ' column contains 1 (Since segmentation output for model has CSF region'
        ' label as 1.')
    parser.add_argument('dist_thres', help='distance from the tissue up to '
                        'which  to fill CSF')
    parser.add_argument(
        'pad_rad', help='Pad radius to use when cropping using the boun'
        'ding box of the region that included newly added CSF regions.')
    parser.add_argument('res_dir', help='path to save results')
    parser.add_argument(
        'fnames_suffix', help='Suffix to add to all output filenames')
    parser.add_argument(
        'model_seg_suffix', help='Suffix to add to the output segmentation file'
        ' created for the model. Useful if you want to know which regions are '
        'considered as CSF regions which depends on the labels present in the '
        'file you provided with -fs_csf_labels. ')
    parser.add_argument(
        'in_new_cluster', action='store_true', help='Env variables for new nef')
    ops = parser.parse_args()
    if not op.exists(ops.res_dir):
        os.makedirs(ops.res_dir)
    return ops

def set_binaries_and_paths(ops):
    """ Set all the required executables in this module and make them global.
    envirionment variable. The paths set for the following environment variables
    are returned (respecting the order) :
    ADLEM_DIR
    ANTS_BIN
    """
    global adlem_dir, ants_dir
    global img_math, binarize, get_dist, img_from_label_img
    adlem_dir = os.getenv('ADLEM_DIR')
    ants_dir = os.getenv('ANTS_BIN')
    if adlem_dir is None:
        raise ValueError('environment variable ADLEM_DIR not found')
    if ants_dir is None:
        raise ValueError('environment variable ANTS_BIN not found.')
    img_math = ants_dir + '/ImageMath'
    if ops.in_new_cluster:
        build_dir = 'buildNewNef'
    else:
        build_dir = 'build'
    binarize = op.join(adlem_dir, build_dir, 'src/BinarizeThreshold')
    get_dist = op.join(adlem_dir, build_dir, 'src/signedDanielssonDistance')
    img_from_label_img = op.join(
        adlem_dir, build_dir, 'src/createImageFromLabelImage')
    return

def crop_by_mask(dim, in_img, out_img, label_mask_img, label='1', padRadius='0'):
    '''
    Run Ants ExtractRegionFromImageByMask command
    '''
    crop = '%s/ExtractRegionFromImageByMask %s' % (ants_dir, dim)
    cmd = ('%s %s %s %s %s %s'
           % (crop, in_img, out_img, label_mask_img, label, padRadius))
    bu.print_and_execute(cmd)
    return


def create_csf_mask_images(ops):
    '''
    Create binary images (masks) with various csf images.
    CSF here can be of two types:
    1. ROIs in ops.in_seg with labels present in ops.fs_csf_labels
    file. Denoted by: FS_CSF
    2. Sulcal CSF computed here using distance map. Denoted by: sulcal_CSF

    Creates image files (that are also returned) are (in the order)
    fs_brain_csf_mask: binary mask with foreground = brain + FS_CSF
    fs_brain_no_csf_mask: binary mask, foreground = brain only
    fs_csf_mask: binary mask, foreground = FS_CSF
    sulcal_csf: binary mask, foreground = sulcal_CSF
    '''
    fl_pref = op.join(ops.res_dir, 'create_csf_mask_images_')
    fs_brain_csf_mask = fl_pref + 'tmp++InBinaryWithCsf.nii.gz'
    cmd = ('%s -i %s -o %s -l 1 -u 15000'
           % (binarize, ops.in_seg, fs_brain_csf_mask))
    bu.print_and_execute(cmd)

    fs_csf_mask = fl_pref + 'tmp++fs_csf_mask.nii.gz'
    cmd = ('%s -t %s -l %s -o %s'
           % (img_from_label_img, ops.fs_csf_labels, ops.in_seg, fs_csf_mask))
    bu.print_and_execute(cmd)
    # Make FS_CSF regions to have label values > largest FS CSF label
    fs_csf_mask_mult = fl_pref + 'tmp++fs_csf_mask_mult.nii.gz'
    # max CSF label = 221. So multiply with say 300 (but -ve)
    cmd = '%s 3 %s m %s -300' % (img_math, fs_csf_mask_mult, fs_csf_mask)
    bu.print_and_execute(cmd)

    # Mask with foreground = only FS segmented tissue and no CSF at all
    fs_brain_no_csf_mask = fl_pref + 'tmp++InBinary.nii.gz'
    cmd = ('%s 3 %s + %s %s'
           % (img_math, fs_brain_no_csf_mask, ops.in_seg, fs_csf_mask_mult))
    bu.print_and_execute(cmd) #Results in -ve values for all CSF regions.
    # Binarize by excluding FS CSF
    cmd = ('%s -i %s -o %s -l 1 -u 15000'
           % (binarize, fs_brain_no_csf_mask, fs_brain_no_csf_mask))
    bu.print_and_execute(cmd)

    # get the distance of the binary input
    dist_img = fl_pref + 'tmp++dist_img.nii.gz'
    cmd = '%s %s %s' % (get_dist, fs_brain_no_csf_mask, dist_img)
    bu.print_and_execute(cmd)
    # Extract desired CSF region input distance threshold
    sulcal_csf_mask = fl_pref + 'tmp++sulcal_csf_mask.nii.gz'
    cmd = ('%s -i %s -o %s -l 0 -u %s'
           % (binarize, dist_img, sulcal_csf_mask, ops.dist_thres))
    bu.print_and_execute(cmd)
    return fs_brain_no_csf_mask, fs_brain_csf_mask, fs_csf_mask, sulcal_csf_mask

def create_seg_for_model(fs_tissue_mask, fs_csf_mask, sulcal_csf_mask, out_seg):
    '''
    Create segmentation image with the following labels:
    Non brain region (NBR) - 0
    CSF - 1  and Tissue - 2
    '''
    # Make tissue as 2
    cmd = '%s 3 %s m %s 2' % (img_math, out_seg, fs_tissue_mask)
    bu.print_and_execute(cmd)
    #if fs_csf_mask is not None: #Include freesurfer csf regions
    cmd = ('%s 3 %s + %s %s'
           % (img_math, out_seg, out_seg, fs_csf_mask))
    bu.print_and_execute(cmd)
    # Now combine sulcal CSF regions with CSF label as 1
    # Add only at those places where freesurfer had zero labels.
    cmd = ('%s 3 %s addtozero %s %s'
           % (img_math, out_seg, out_seg, sulcal_csf_mask))
    bu.print_and_execute(cmd)


def create_seg_for_atrophy(ops, sulcal_csf_mask, out_seg):
    '''
    Create segmentation image with all the FreeSurfer labels intact +
    additional sulcal CSF added where the label asssinged to the CSF
    is the same as the CSF label of FreeSurfer (24)
    '''
    fl_pref = op.join(ops.res_dir, 'create_seg_for_atrophy_')
    # Make sulcal CSF regions to have label 24.
    sulcal_csf_mask_mult = fl_pref + 'tmp++sulcal_csf_mask_mult.nii.gz'
    cmd = '%s 3 %s m %s 24' % (img_math, sulcal_csf_mask_mult, sulcal_csf_mask)
    bu.print_and_execute(cmd)
    # Add this to input FS seg only at those places where it had zero labels.
    cmd = ('%s 3 %s addtozero %s %s'
           % (img_math, out_seg, ops.in_seg, sulcal_csf_mask_mult))
    bu.print_and_execute(cmd)


def create_crop_mask(fs_brain_csf_mask, sulcal_csf_mask, out_seg):
    '''
    Create output crop mask.
    Mask foreground = sulcal CSF regions + foreground of fs_brain_csf_mask.
    '''
    # Add sulcal_csf_mask to fs_brain_csf_mask only where fs_brain_csf_mask = 0
    cmd = ('%s 3 %s addtozero %s %s'
           % (img_math, out_seg, fs_brain_csf_mask, sulcal_csf_mask))
    bu.print_and_execute(cmd)


def main():
    """ Generate segmented images, mask and atrophy map that are adapted to be
    inputs for AdLemModel.
    """
    ops = get_input_options()
    set_binaries_and_paths(ops)
    # Output filenames
    crop_info = 'D%sR%s' % (ops.dist_thres, ops.pad_rad)
    out_seg_for_model = op.join(
        ops.res_dir, 'segForModel%s%s%s.nii.gz'
        % (crop_info, ops.fnames_suffix, ops.model_seg_suffix))
    out_seg_for_atrophy = op.join(
        ops.res_dir, 'labelsWithCsfForAtrophy%s%s.nii.gz'
        % (crop_info, ops.fnames_suffix))
    out_t1 = op.join(ops.res_dir,
                     't1%s%s.nii.gz' % (crop_info, ops.fnames_suffix))
    out_crop_mask = op.join(ops.res_dir, 'cropMaskD%s%s.nii.gz'
                            % (ops.dist_thres, ops.fnames_suffix))

    # Create different temporary mask images
    (fs_brain_no_csf_mask, fs_brain_csf_mask,
     fs_csf_mask, sulcal_csf_mask) = create_csf_mask_images(ops)

    create_seg_for_model(fs_brain_no_csf_mask, fs_csf_mask, sulcal_csf_mask,
                         out_seg_for_model)
    # Create output segmentation image for atrophy generation purposes.
    create_seg_for_atrophy(ops, sulcal_csf_mask, out_seg_for_atrophy)
    # create crop mask
    create_crop_mask(fs_brain_csf_mask, sulcal_csf_mask, out_crop_mask)

    # Crop all the images:
    crop_by_mask(3, out_seg_for_model, out_seg_for_model, out_crop_mask,
                 padRadius=ops.pad_rad)
    crop_by_mask(3, out_seg_for_atrophy, out_seg_for_atrophy,
                 out_crop_mask, padRadius=ops.pad_rad)
    crop_by_mask(3, ops.t1_img, out_t1, out_crop_mask, padRadius=ops.pad_rad)

    # Delete all tmp++ prefixed files from the result directory.
    cmd = 'rm ' + op.join(ops.res_dir, '*tmp++*')
    bu.print_and_execute(cmd)


if __name__ == "__main__":
    main()
