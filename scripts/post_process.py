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
    parser.add_argument(
        'config_file', help='input config file which contains options database '
        'as a python dictionary.')
    parser.add_argument(
        '--paste_image', action='store_true', help='Paste the cropped images '
        'back to original size')
    parser.add_argument(
        '--resample_image', action='store_true', help='Resample intensity from '
        'another image for the simulated image.')
    parser.add_argument(
        '-s', '--sim_prefixes', help=' Sim prefixes used in simulating images. '
        'Separate each of them with comma. If not given, searches it in config '
        'where it must be provided as a list. If given, sim_prefixes key in the '
        'config file is ignored.')
    parser.add_argument(
        '-p', '--patients', help=' Patient ids. Separate them with comma. If '
        'not given, searches it in config where it must be provided as a list. '
        ' If given, sim_prefixes key in the config file is ignored.')
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


def launch_job(job_string, job_name, res_dir, is_new_nef = False, print_only = False):
    """ Launch post processing from config file.
    """
    cmd = bu.sophia_nef_pbs_setting(is_new_nef) + job_string
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
        adlem_dir, 'buildNewNef/src/pasteImageToBiggerImage')
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

def paste_image(ops_dict, ops, patient):
    '''
    Paste cropped image to bigger size using paste_img binary.
    '''
    global cmd_all
    d, p = ops_dict, patient
    big_img_dir = op.join(d['data_base_dir'], p, d['in_img_dir'])
    big_img_suffix = p[d['in_img_pid_skip']:] + d['in_img_suffix']
    big_img = op.join(big_img_dir, d['in_img_prefix']+big_img_suffix)
    sim_dir = op.join(d['sim_base_dir'], p)
    crop_mask = op.join(sim_dir, d['sim_crop_mask'])
    for sim_prefix in d['sim_prefixes']:
        in_imgs = [op.join(sim_dir, sim_prefix+x)
                   for x in d['sim_cropped_img_suffices']]
        out_imgs = [op.join(sim_dir, sim_prefix+x)
                    for x in d['sim_full_img_suffices']]
        for in_img, out_img in zip(in_imgs, out_imgs):
            cmd = ('%s %s %s %s %s\n\n'
                   % (paste_img, in_img, big_img, crop_mask, out_img))
            #print(cmd)
            if op.exists(in_img) and (not op.exists(out_img) or ops.overwrite):
                cmd_all = bu.update_or_execute_cmd(
                    cmd_all, cmd, not ops.in_cluster)
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
    data_dir = op.join(d['data_base_dir'], p)
    in_img_dir = op.join(data_dir, d['in_img_dir'])
    pid = p[d['in_img_pid_skip']:]
    rl_img_suffix = pid + d['in_img_suffix']
    ref_img = op.join(in_img_dir, d['in_img_prefix']+rl_img_suffix)
    mvng_img = op.join(in_img_dir, d['in_resample_img_prefix']+rl_img_suffix)
    reg_disp_field = op.join(data_dir, d['reg_dir'],
                             d['disp_field_prefix']+pid+d['svf_suffix'])
    sim_dir = op.join(d['sim_base_dir'], p)
    for sim_prefix in d['sim_prefixes']:
        sim_vels = [op.join(sim_dir, sim_prefix+x)
                    for x in d['sim_vel_suffices']]
        out_imgs = [op.join(sim_dir, sim_prefix+x)
                    for x in d['sim_resampled_img_suffices']]
        warp_fields = [op.join(sim_dir, sim_prefix+x)
                       for x in d['resampling_disp_suffices']]
        for sim_vel, out_img, warp_field in zip(
                sim_vels, out_imgs, warp_fields):
            if op.exists(sim_vel) and (not op.exists(out_img) or ops.overwrite):
                cmd = ('%s 3 %s -R %s %s %s\n\n'
                       % (compose_transform, warp_field,
                          ref_img, sim_vel, reg_disp_field))
                cmd_all = bu.update_or_execute_cmd(
                    cmd_all, cmd, not ops.in_cluster)
                cmd = ('%s -d 3 -e 0 -i %s -r %s -o %s -n BSpline -t %s\n\n' %
                       (apply_transform, mvng_img, ref_img,
                        out_img, warp_field))
                cmd_all = bu.update_or_execute_cmd(
                    cmd_all, cmd, not ops.in_cluster)
                #print(cmd)

def main():
    """ Post process images . Run with -h to see options."""
    global cmd_all
    ops = get_input_options()
    set_binaries_and_paths()
    cmd_all = vars_setup_string()
    d = bu.read_dict_from_file(ops.config_file)
    if ops.patients:
        d['patients'] = ops.patients.split(',')
    if ops.sim_prefixes:
        d['sim_prefixes'] = ops.sim_prefixes.split(',')
    #print d
    for patient in d['patients']:
        if ops.paste_image:
            paste_image(d, ops, patient)
        if ops.resample_image:
            resample_intensity(d, ops, patient)
    if ops.in_cluster:
        if not ops.job_name or not ops.log_dir:
            raise ValueError('Provide both -l[--log_dir] and -j[--job_name] '
                             'when using --in_cluster')
        launch_job(cmd_all, ops.job_name, ops.log_dir)
        #launch_job(cmd_all, ops.job_name, ops.log_dir, False, True)

if __name__ == "__main__":
    main()
