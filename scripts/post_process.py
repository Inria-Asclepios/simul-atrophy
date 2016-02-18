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

# Paths of required binaries
adlem_dir, ants_dir = None, None
# Required binaries
paste_img, compose_transform, apply_transform = None, None, None
#img_math, binarize = None, None
cmd_all = None

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
    parser.add_argument('--in_cluster', action='store_true', help='Run in '
                        'cluster if the option is provided')
    parser.add_argument('-l', '--log_dir', help='Job log dir.')
    parser.add_argument('-j', '--job_name', help='If --in_cluster, provide a '
                        'name of the job to be launched.')
    ops = parser.parse_args()
    return ops

def vars_setup_string():
    """ Returns a set of commands to set required env variables.
    """
    setup_string = """
export ADLEM_DIR=%s
export ANTS_BIN=%s
export PATH=$PATH:$ANTS_BIN
""" % (adlem_dir, ants_dir)
    return setup_string


def launch_job(job_string, job_name, res_dir, print_only = False):
    """ Launch post processing from config file.
    """
    cmd = bu.sophia_nef_pbs_setting() + job_string
    if print_only:
        print(cmd)
    else:
        bu.qsub_job(name=job_name, queue='seqshort', walltime='03:00:00',
                    procs='nodes=1:xeon:ppn=1', mem='mem=8gb',
                    dest_dir=res_dir, add_atrb=None, cmd=cmd)
    return

def set_binaries_and_paths():
    """ Set all the required executables in this module and make them global.
    envirionment variable. The paths set for the following environment variables
    are returned (respecting the order) :
    ADLEM_DIR
    ANTS_BIN
    """
    global adlem_dir, ants_dir
    global paste_img, compose_transform, apply_transform
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
    compose_transform = op.join(ants_dir, 'ComposeMultiTransform')
    apply_transform = op.join(ants_dir, 'antsApplyTransforms')
    return

def svf_exp(euler_integration, scale, in_svf, out_disp_field, in_cluster):
    '''
    Interface to svf exponential integration. Returns command line string.
    Computes and saves out_disp_field = exp(scale X in_svf)
    <euler_integration>: 1=>Euler itegration; 0=>scaling and squaring
    '''
    if in_cluster is None:
        binary = op.join('/home/bkhanal/works/tools/svfExp', 'SVFExponential')
    else:
        binary = op.join('/epi/asclepios2/bkhanal/works/tools/svfExp',
                         'SVFExponential')
    cmd = ('%s -e %s -s %s -o %s -i %s\n\n' %
           (binary, euler_integration, scale, out_disp_field, in_svf))
    return cmd

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
    global cmd_all
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
        if op.exists(in_img) and (not op.exists(out_img) or ops.overwrite):
            crop_mask = op.join(sim_dir, d['crop_img'])
            pid = p.partition('_')[2]
            if d['seg_file_of'] == 'baseline':
                big_img = op.join(
                    data_dir, d['real_img_dir'], d['bl_prefix']+pid+ext)
            else:
                big_img = op.join(
                    data_dir, d['real_img_dir'], d['fl_prefix']+pid+ext)
            cmd = ('%s %s %s %s %s\n\n'
                   % (paste_img, in_img, big_img, crop_mask, out_img))
            cmd_all = bu.update_or_execute_cmd(cmd_all, cmd, not ops.in_cluster)
            #print(cmd)

def resample_intensity(ops_dict, ops, patient):
    '''
    Resample Intensity from the desired image (different than the baseline used
    to obtain the simulation transormation field sim_vel.). This is done by
    composing the transformation field as follows:
       M[T_r(T_nr(T_s_inv(x)))]   where,
    M is the image from which to resample intensity
    T_r Transformation obtained from rigid registration between M and B where
    B is the input image that was used to simulate the transformation field T_s.
    T_nr Transformation obtained from non-rigid registration between B and
    M aligned to B wth T_r.
    If M and B were already rigidly aligned, T_r is identity.
    '''
    global cmd_all
    d, p = ops_dict, patient
    pid = p.partition('_')[2]
    data_dir = op.join(d['data_base_dir'], patient)
    sim_dir = op.join(d['sim_base_dir'], patient)
    sim_prefix, ext = ops.sim_prefix, '.nii.gz'
    tpts = d['tpts'].split(',')
    for t in tpts:
        if d['seg_file_of'] == 'baseline':
            ref_img = op.join(
                data_dir, d['real_img_dir'], d['bl_prefix']+pid+ext)
            mvng_img = op.join(
                data_dir, d['real_img_dir'], d['fl_prefix']+pid+ext)
            reg_disp_field = op.join(data_dir, d['reg_dir'],
                                     d['disp_field_prefix']+pid+d['svf_suffix'])
        else:
            ref_img = op.join(
                data_dir, d['real_img_dir'], d['fl_prefix']+pid+ext)
            mvng_img = op.join(
                data_dir, d['real_img_dir'], d['bl_prefix']+pid+ext)
            reg_disp_field = op.join(
                data_dir, d['reg_dir'],
                d['minus_disp_field_prefix']+pid+d['svf_suffix'])
        sim_vel = op.join(sim_dir, '%sT%svel%s' % (sim_prefix, t, ext))
        out_img = op.join(sim_dir, '%sFullSizeSimulatedResampledT%s%s'
                          % (sim_prefix, t, ext))
        if op.exists(sim_vel) and (not op.exists(out_img) or ops.overwrite):
            warp_field = op.join(sim_dir, '%sFullSizeResampleDispFieldT%s%s'
                                 % (sim_prefix, t, ext))
            cmd = ('%s 3 %s -R %s %s %s\n\n'
                   % (compose_transform, warp_field,
                      ref_img, sim_vel, reg_disp_field))
            cmd_all = bu.update_or_execute_cmd(cmd_all, cmd, not ops.in_cluster)
            # Resample intensity from the image not used for simulation.
            cmd = ('%s -d 3 -e 0 -i %s -r %s -o %s -n BSpline -t %s\n\n' %
                   (apply_transform, mvng_img, ref_img, out_img, warp_field))
            cmd_all = bu.update_or_execute_cmd(cmd_all, cmd, not ops.in_cluster)
            #print(cmd)


def main():
    """ Generate atrophy maps. Run with -h to see options."""
    global cmd_all
    ops = get_input_options()
    set_binaries_and_paths()
    cmd_all = vars_setup_string()
    d = parse_config_file(ops.config_file) # dictionary with options database
    d['patients'] = d['patients'].split(',')
    #print d
    for patient in d['patients']:
        paste_image(d, ops, patient)
        resample_intensity(d, ops, patient)
    if ops.in_cluster:
        if not ops.job_name or not ops.log_dir:
            raise ValueError('Provide both -l[--log_dir] and -j[--job_name] '
                             'when using --in_cluster')
        launch_job(cmd_all, ops.job_name, ops.log_dir)
        #launch_job(cmd_all, ops.job_name, ops.log_dir, True)

if __name__ == "__main__":
    main()
