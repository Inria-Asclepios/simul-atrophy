#!/usr/bin/env python
import subprocess
import argparse as ag
import os
import os.path as op
import time as tm
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
    parser.add_argument('in_seg', help='input freesurfer aparc+aseg file in '
                        'the same space as the input MR image.')
    parser.add_argument('dist_thres', help='distance from the tissue up to '
                        'which  to fill CSF')
    parser.add_argument(
        'pad_rad', help='Pad radius to use when cropping using the boun'
        'ding box of the region that included newly added CSF regions.')
    parser.add_argument(
        'atrophy_tables', help='filenames of the tables having label values and'
        'the corresponding desired atrophy values. Separate individual files '
        'with comma WITHOUT space.')
    parser.add_argument('res_dir', help='path to save results')
    parser.add_argument(
        '--only_sulcal_csf', action='store_true', help='If given considers '
        'all the ventricles as tissue + all CSF segmented by FS as tissue. Only'
        ' sulcal CSF computed in this script will be considered as CSF and the '
        ' output segmentation for the model will be created accordingly. That '
        'is all the ventricles will be set the same label as that of tissue')
    ops = parser.parse_args()
    if not op.exists(ops.res_dir):
        os.makedirs(ops.res_dir)
    return ops

def set_binaries_and_paths():
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
    binarize = adlem_dir + '/build/src/BinarizeThreshold'
    get_dist = adlem_dir + '/build/src/signedDanielssonDistance'
    img_from_label_img = adlem_dir + '/build/src/createImageFromLabelImage'
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
    Create image files (that are also returned) are (in the order)
    fs_brain_csf_mask: binary mask with foreground = brain + csf
    fs_brain_no_csf_mask: binary mask, foreground = brain only
    fs_csf_mask: binary mask, foreground = FS segmented csf only
    sulcal_csf: binary mask, foreground = sulcal CSF computed using distance map
    '''
    # Binarize input FreeSurfer segmentation in following two ways
    # 1. with brain + FS segmeted CSF : fs_brain_csf_mask
    # 2. Only tissue excluding CSF (fs_brain_no_csf_mask)
    # The sulcal CSF region is obtained by using distance map on fs_brain_no_csf_mask
    # Mask with foreground = all FS brain and CSF segmented regions.
    fl_pref = op.join(ops.res_dir, 'create_csf_mask_images_')
    fs_brain_csf_mask = fl_pref + 'tmp++InBinaryWithCsf.nii.gz'
    cmd = ('%s -i %s -o %s -l 1 -u 15000'
           % (binarize, ops.in_seg, fs_brain_csf_mask))
    bu.print_and_execute(cmd)

    # Mask with foreground = Only FS segmented CSF regions.
    if ops.only_sulcal_csf:
        fs_csf_mask = None
        fs_brain_no_csf_mask = None
    else:
        fs_csf_mask = fl_pref + 'tmp++fs_csf_mask.nii.gz'
        cmd = ('%s -t %s/configFiles/freeSurferCsfLabels -l %s -o %s'
               % (img_from_label_img, adlem_dir, ops.in_seg, fs_csf_mask))
        bu.print_and_execute(cmd)
        # Make FS CSF regions to have label values > largest FS CSF label
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
    if ops.only_sulcal_csf: # Here we consider all FS segmented CSF as tissue
        cmd = '%s %s %s' % (get_dist, fs_brain_csf_mask, dist_img)
    else: # Here FS segmented CSF (and ventricles are taken as CSF)
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
    if fs_csf_mask is not None: #Include freesurfer csf regions
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


def create_atrophy_map(ops, seg_for_atrophy, out_atrophy_map):
    '''
    Create atrophy map from input tables and the input segmentation image.
    '''
    # Now create atrophy map: The first atrophy table is used to create while
    # the subsequent ones modify the one created with the first one.
    cmd_pref = ('%s -l %s -o %s -t '
                % (img_from_label_img, seg_for_atrophy, out_atrophy_map))
    tables = ops.atrophy_tables.split(',') # files separated by comma
    # print('list of tables: ')
    # print tables
    cmd = cmd_pref + tables[0]
    bu.print_and_execute(cmd)
    for table in tables[1:]: # from the remaining tables modify the output file
        cmd = cmd_pref + table + ' -m ' + out_atrophy_map
        bu.print_and_execute(cmd)
        tm.sleep(0.2) # wait cmd to write file before launching another one.


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
    set_binaries_and_paths()
    # Output filenames
    crop_info = 'D%sR%s' % (ops.dist_thres, ops.pad_rad)
    out_seg_for_model = op.join(ops.res_dir, 'maskwithCsf%s.nii.gz'
                                % (crop_info))
    out_seg_for_atrophy = op.join(ops.res_dir, 'labelsWithCsf%s.nii.gz'
                                  % (crop_info))
    out_atrophy = op.join(ops.res_dir, 'atrophyMap%s.nii.gz' % (crop_info))
    out_t1 = op.join(ops.res_dir, 't1%s.nii.gz' % (crop_info))
    out_crop_mask = op.join(ops.res_dir, 'cropMaskD%sPadRadius%s.nii.gz'
                            % (ops.dist_thres, ops.pad_rad))

    # Create different temporary mask images
    (fs_brain_no_csf_mask, fs_brain_csf_mask,
     fs_csf_mask, sulcal_csf_mask) = create_csf_mask_images(ops)

    # Create output segmentation image for the model
    if ops.only_sulcal_csf:
        fs_tissue_mask = fs_brain_csf_mask
    else:
        fs_tissue_mask = fs_brain_no_csf_mask
    create_seg_for_model(fs_tissue_mask, fs_csf_mask, sulcal_csf_mask,
                         out_seg_for_model)
    # Create output segmentation image for atrophy generation purposes.
    create_seg_for_atrophy(ops, sulcal_csf_mask, out_seg_for_atrophy)
    # Create atrophy map:
    create_atrophy_map(ops, out_seg_for_atrophy, out_atrophy)
    # create crop mask
    create_crop_mask(fs_brain_csf_mask, sulcal_csf_mask, out_crop_mask)

    # Crop all the images:
    crop_by_mask(3, out_atrophy, out_atrophy, out_crop_mask, padRadius='0')
    crop_by_mask(3, out_seg_for_model, out_seg_for_model, out_crop_mask,
                 padRadius='0')
    crop_by_mask(3, out_seg_for_atrophy, out_seg_for_atrophy,
                 out_crop_mask, padRadius='0')
    crop_by_mask(3, ops.t1_img, out_t1, out_crop_mask, padRadius='0')

    # Delete all tmp++ prefixed files from the result directory.
    cmd = 'rm ' + op.join(ops.res_dir, '*tmp++*')
    bu.print_and_execute(cmd)


if __name__ == "__main__":
    main()
