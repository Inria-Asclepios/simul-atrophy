#!/usr/bin/env python
import subprocess
import sys
import argparse as ag
import bish_utils as bu

def get_input_options():
    ''' Command line interface, get user input options.
    '''
    parser = ag.ArgumentParser()
    parser.add_argument('inFile', help='Input Image to be warped')
    parser.add_argument('velFile', help='Input vel-field')
    parser.add_argument('outPref', help='warped output filenames prefix')
    parser.add_argument('invertField', help='true or false: invert the field '
                        'before warping or not')
    parser.add_argument('numOfScales', help='Total Number of time points to '
                        ' warp image to', type=int)
    parser.add_argument('firstScale', help='Starting time point', type=int)
    parser.add_argument('stepSize', help='Step size between two time points',
                        type=float)
    parser.add_argument(
        'compScheme', help='Computational scheme: 0 for scaling and squaring, '
        'is faster; 1 for forward euler, is more accurate.')
    parser.add_argument('interp', help='Possible types: '
                        'linear/bspline/nearestneighbor')
    parser.add_argument('-bspline_order', help='When bspline is used')
    parser.add_argument(
        '--save_only_4D', action='store_true', help='save only 4D image.')

    ops = parser.parse_args()
    print('%s will be warped\n' % (ops.inFile,))
    print('number of timesteps: %d\n' % (ops.numOfScales,))
    print('first scale value: %d\n' % (ops.firstScale,))
    print ('step size: %f\n' % (ops.stepSize,))
    return ops

def main():
    '''
    Warp input image using the input velocity field at different scales. The
    scales are integrated by considering input field as a SVF and using either
    forward Euler or scaling and squaring.
    '''
    #executables
    work_dir = "/user/bkhanal/home/works/"
    svf_exp = work_dir + "tools/svfExponential/SVFExponential"
    warp_img = work_dir + "AdLemModel/build/src/WarpImage"
    combine_imgs = work_dir + "tools/combineImages/build/combineImages"
    ext = '.nii.gz'

    #inputs
    ops = get_input_options()
    in_img = ops.inFile
    in_vel = ops.velFile
    scale = ops.firstScale
    num = 1
    out_img_pref = ops.outPref+'WarpedScale'
    out_4D = '%sWarped4D%sSteps%s' % (ops.outPref, ops.numOfScales, ext)
    # Make sure that out_img_pref string is not the starting string of out_4D.
    # Because I delete individual files with out_img_pref*ext !!
    if ops.interp == 'bspline':
        interp_ops = ('--interpolator %s --order %s' % (
            ops.interp, ops.bspline_order))
    else:
        interp_ops = '--interpolator %s' % (ops.interp)
    while num <= ops.numOfScales:
        out_vel = '%svelScaled%s%s' % (ops.outPref, scale, ext)
        out_img = '%s%s%s' % (out_img_pref, num, ext)
        #create vel field from in_vel
        cmd = ('%s -s %s -e %s -i %s -o %s' % (
            svf_exp, scale, ops.compScheme, in_vel, out_vel))
        bu.print_and_execute(cmd)
        #create warped Image from in_img using out_vel created with svf_exp
        cmd = ('%s --inImage %s --displacementImage %s '
               '--outImage %s --invert %s %s' % (
                   warp_img, in_img, out_vel, out_img, ops.invertField,
                   interp_ops))
        bu.print_and_execute(cmd)
        #Delete vel field
        bu.print_and_execute('rm ' + out_vel + '\n')
        #Go to next file
        num += 1
        scale += ops.stepSize
    #combine images: ops: numOfScales imageType imageBaseFileName out_imgFileName
    cmd = ('%s %s %s %s %s' % (
        combine_imgs, ops.numOfScales, ext, out_img_pref, out_4D))
    bu.print_and_execute(cmd)
    #delete individual warped images:
    if ops.save_only_4D:
        bu.print_and_execute('rm %s*%s' % (out_img_pref, ext))

if __name__ == "__main__":
    main()
