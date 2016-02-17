#!/usr/bin/python
""" Post process after simulating the images.
"""
# Author: Bishesh Khanal <bishesh.khanal@inria.fr>
#
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
paste_img = None
#img_math, binarize = None, None


def get_input_options():
    """ Command line interface, get input options and interact with the user.
    """
    parser = ag.ArgumentParser('Post process the output simulated results.'
                               ' using given config file.')
    parser.add_argument('config_file', help='input config file')
    parser.add_argument('sim_prefix', help=' prefix used in simulating images')
    parser.add_argument(
        '--overwrite', action='store_true', help='Overwrite the results even if'
        ' existing files are present. ')
    ops = parser.parse_args()
    return ops

def set_binaries_and_paths():
    """ Set all the required executables in this module and make them global.
    envirionment variable. The paths set for the following environment variables
    are returned (respecting the order) :
    ADLEM_DIR
    ANTS_BIN
    """
    global adlem_dir, ants_dir
    global paste_img
    adlem_dir = os.getenv('ADLEM_DIR')
    ants_dir = os.getenv('ANTS_BIN')
    if adlem_dir is None:
        raise ValueError('environment variable ADLEM_DIR not found')
    if ants_dir is None:
        raise ValueError('environment variable ANTS_BIN not found.')
    #img_math = op.join(ants_dir, 'ImageMath')
    #binarize = op.join(adlem_dir, 'build/src/BinarizeThreshold')
    paste_img = op.join(
        adlem_dir, 'build/src/pasteImageToBiggerImage')
    return


def parse_config_file(config_file):
    '''
    Parse the config file and return a dictionary.
    '''
    cfg_keys, cfg_vals = [], []
    with open(config_file, 'r') as fl:
        for line in fl:
            if line.startswith('#'):
                cfg_keys += [x for x in line.split() if '#' not in x]
            else:
                cfg_vals += line.split()
        #cfg_keys = [line.strip() for line in fl if line.startswith('#')]
    cfg_dict = dict(zip(cfg_keys, cfg_vals))
    #cfg_dict['patients'] = cfg_dict['patients'].split(',')
    return cfg_dict

def paste_image(ops_dict, ops, patient):
    '''
    Paste cropped image to bigger size using paste_img binary.
    '''
    d, p = ops_dict, patient
    data_dir = op.join(d['data_base_dir'], patient)
    sim_dir = op.join(d['sim_base_dir'], patient)
    sim_prefix, ext = ops.sim_prefix, '.nii.gz'
    tpts = d['tpts'].split(',')
    for t in tpts:
        out_img = op.join(sim_dir, '%sFullSizeSimulatedT%s%s'
                                  % (sim_prefix, t, ext))
        in_img = op.join(sim_dir, '%sWarpedImageBsplineT%s%s'
                         % (sim_prefix, t, ext))
        if op.exists(in_img) and not op.exists(out_img):
            crop_mask = op.join(sim_dir, d['crop_img'])
            pid = p.partition('_')[2]
            if d['seg_file_of'] == 'baseline':
                big_img = op.join(
                    data_dir, d['real_img_dir'], d['bl_prefix']+pid+ext)
            else:
                big_img = op.join(
                    data_dir, d['real_img_dir'], d['fl_prefix']+pid+ext)
            cmd = ('%s %s %s %s %s' % (paste_img, in_img, big_img,
                                       crop_mask, out_img))
            bu.print_and_execute(cmd)
            #print(cmd)


def main():
    """ Generate atrophy maps. Run with -h to see options."""
    ops = get_input_options()
    set_binaries_and_paths()
    d = parse_config_file(ops.config_file) # dictionary with options database
    d['patients'] = d['patients'].split(',')
    #print d
    data_dir, sim_dir = d['data_base_dir'], d['sim_base_dir']
    sim_prefix, ext = ops.sim_prefix, '.nii.gz'
    for patient in d['patients']:
        paste_image(d, ops, patient)


if __name__ == "__main__":
    main()
