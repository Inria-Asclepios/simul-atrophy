#!/usr/bin/python
""" Create atrophy map from segmentation label images, and input atrophy tables
or input images.
"""
# Author: Bishesh Khanal <bishesh.khanal@inria.fr>
# Asclepios INRIA Sophia Antipolis
#

import argparse as ag
import os
import os.path as op
import time as tm
import bish_utils as bu
import image_math as im

# Paths of required binaries
adlem_dir, ants_dir = None, None
# Required binaries
uniformize_regionally = None
#img_math, binarize = None, None
img_from_label_img = None


def get_input_options():
    """ Command line interface, get input options and interact with the user.
    """
    parser = ag.ArgumentParser('Write atrophy map either from given tables'
                               ' or from input intensity image.')
    parser.add_argument('out_atrophy', help='output atrophy map image')
    parser.add_argument(
        '-seg_for_atrophy', help='Input label image with all ROIs where atrophy '
        'values are to be set. If you provide -atrophy_img and -crop_img, Make '
        'sure that seg_for_atrophy has the same size as that of the one output '
        'atrophy map will have when cropped with given options')
    parser.add_argument(
        '-atrophy_tables', help='filenames of the tables having label values '
        'and the corresponding desired atrophy values. Separate individual '
        'files with comma WITHOUT space. Tables from this option are used at '
        'very end. So for e.g. if you use -atrophy_img also, the values of '
        'output atrophy computed can be modified by this table!')
    parser.add_argument(
        '-atrophy_img', help='Input image file to adapt as atrophy map. If not'
        ' given, starts from image with 0 everywhere. The tables with option '
        '-atrophy_tables will be used after using the atrophy_img and related '
        'other options.')
    parser.add_argument(
        '--negate_atrophy_vals', action='store_true', help='Input atrophy vals '
        'from  will be multiplied by -1 for the final output.')
    atrophy_img = parser.add_argument_group('When using -atrophy_img')
    atrophy_img.add_argument('-crop_mask', help='If given, I will crop input '
                             '-atrophy_img before using it.')
    atrophy_img.add_argument('-crop_pad_rad', help='Pad radius to use when '
                             'cropping input -atrophy_img.')
    atrophy_img.add_argument(
        '--uniform_regional_means', action='store_true', help='Compute and set '
        'uniform mean intensity of -atrophy_img in each ROI. ROIs are obtained '
        'from the -uniformize_labels and seg_for_atrophy.')
    atrophy_img.add_argument(
        '-uniformize_labels', help='File with all the labels defining the ROIs '
        'in seg_for_atrophy to compute mean. Relevant only when --uniform_'
        'regional_means is provided.')
    atrophy_img.add_argument(
        '-seg_for_uniformize', help='Input label image with all the labels '
        'defining the ROIs present in -uniformize_labels to compute mean. '
        'Relevant only when --uniform_regional_means is provided.')
    atrophy_img.add_argument(
        '--start_from_atrophy_img', action='store_true', help='If given, '
        'uniformized image will have  values copied from atrophy_img for all '
        'ROIs not covered by -uniformize labels; otherwise these non selected '
        'regions will be set to zero value.')
    ops = parser.parse_args()
    if (ops.atrophy_tables and not ops.seg_for_atrophy) or (ops.seg_for_atrophy and not ops.atrophy_tables):
        parser.error('Either both or none of -seg_for_atrophy and '
                     '-atrophy_tables options are required')
    if ops.uniform_regional_means and not ops.uniformize_labels:
        parser.error('-uniformize_labels needed for --uniform_regional_means')
    if ops.uniformize_labels and not ops.uniform_regional_means:
        parser.error('Cannot use -uniformize_labels when '
                     '--uniform_regional_means not provided')
    if ops.atrophy_img and ops.crop_mask:
        if not ops.crop_pad_rad:
            parser.error('Must provide -crop_pad_rad also when -crop_mask is '
                         'provided!')
    return ops

def set_binaries_and_paths():
    """ Set all the required executables in this module and make them global.
    envirionment variable. The paths set for the following environment variables
    are returned (respecting the order) :
    ADLEM_DIR
    ANTS_BIN
    """
    global adlem_dir, ants_dir
    global img_from_label_img, uniformize_regionally #,img_math, binarize
    adlem_dir = os.getenv('ADLEM_DIR')
    ants_dir = os.getenv('ANTS_BIN')
    if adlem_dir is None:
        raise ValueError('environment variable ADLEM_DIR not found')
    if ants_dir is None:
        raise ValueError('environment variable ANTS_BIN not found.')
    #img_math = op.join(ants_dir, 'ImageMath')
    #binarize = op.join(adlem_dir, 'build/src/BinarizeThreshold')
    img_from_label_img = op.join(
        adlem_dir, 'build/src/createImageFromLabelImage')
    uniformize_regionally = op.join(
        adlem_dir, 'build/src/uniformizeWithRegionalMeans')
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

def atrophy_map_from_tables(ops):
    '''
    Create atrophy map from input tables and the input segmentation image.
    Returns created output atrophy image filename.
    '''
    # Now create atrophy map: The first atrophy table is used to create while
    # the subsequent ones modify the one created with the first one.
    cmd_pref = ('%s -l %s -o %s -t '
                % (img_from_label_img, ops.seg_for_atrophy, ops.out_atrophy))
    tables = ops.atrophy_tables.split(',') # files separated by comma
    # print('list of tables: ')
    # print tables
    cmd = cmd_pref + tables[0]
    if ops.atrophy_img:
        cmd += ' -m ' + ops.out_atrophy
    bu.print_and_execute(cmd)
    for table in tables[1:]: # from the remaining tables modify the output file
        cmd = cmd_pref + table
        if not ops.atrohy_img:
            cmd += ' -m ' + ops.out_atrophy
        bu.print_and_execute(cmd)
        tm.sleep(0.2) # wait cmd to write file before writing another one.
    return ops.out_atrophy



def atrophy_map_from_img(ops):
    '''
    Create atrophy map from input image.
    '''
    if not ops.crop_mask and not ops.uniform_regional_means:
        cmd = 'cp %s %s' % (ops.atrophy_img, ops.out_atrophy)
        return
    in_atrophy = ops.atrophy_img
    if ops.crop_mask:
        crop_by_mask(3, in_atrophy, ops.out_atrophy,
                     ops.crop_mask, padRadius=ops.crop_pad_rad)
        in_atrophy = ops.out_atrophy
    if ops.uniform_regional_means:
        cmd = ('%s -t %s -i %s -l %s -u %s'
               % (uniformize_regionally, ops.uniformize_labels, in_atrophy,
                  ops.seg_for_uniformize, ops.out_atrophy))
        if ops.start_from_atrophy_img:
            cmd += ' -s'
        bu.print_and_execute(cmd)
    return

def main():
    """ Generate atrophy maps. Run with -h to see options."""
    ops = get_input_options()
    set_binaries_and_paths()
    if ops.atrophy_img:
        atrophy_map_from_img(ops)
    if ops.atrophy_tables:
        atrophy_map_from_tables(ops)
    if ops.negate_atrophy_vals:
        im.mult_scalar_imgs([ops.out_atrophy], [ops.out_atrophy], -1)

if __name__ == "__main__":
    main()
