#!/usr/bin/env python
import subprocess
import sys
import argparse as ag
import bish_utils as bu
#import time as tm

def get_input_options():
    ''' Command line interface, get user input options.
    '''
    parser = ag.ArgumentParser()
    parser.add_argument('in_img', help='Filename of the input 4D image')
    parser.add_argument('out_gif', help='Filename of the output file')
    parser.add_argument('total_tpts', help='Total number of time points to '
                        'extract from the input image', type=int)
    parser.add_argument('slice_axis', help='Slice axis: 0,1 or 2')
    parser.add_argument(dest='slice_num', help='slice number to be extracted')
    parser.add_argument('resize', help='resize to maximum possible size in each'
                        ' direction, e.g. 400x400')
    parser.add_argument('-time_unit', help='unit of time for display. E.g yrs'
                        '\n If none, time info not overlayed over the video')
    parser.add_argument(
        '-time_step', help='time between tpt1 and tpt2. Used only when '
        'time_unit is also used. Default 1.0', type=float)
    parser.add_argument(
        'delay', help='time delay between frames in milliseconds')
    parser.add_argument(
        '-rotate', help='If given overrides the default rotation used for '
        'proper orientation of the slice.')
    parser.add_argument(
        '-crop', help='convert -crop ops: If given crops 2D slices before '
        'combining. wxh+x+y')

    ops = parser.parse_args()
    return ops


def main():
    '''
    Extract 3D images from a 4D image and then extract one selected slice
    from each of these 3D images and combine them as a gif.
    '''
    ops = get_input_options()
    #split 4d image to get all the 3D images:
    tmp3d_basename, file_ext = 'tmp3DImage', '.nii.gz'
    cmd = 'ImageMath 4 %s%s TimeSeriesDisassemble %s' % (
        tmp3d_basename, file_ext, ops.in_img)
    bu.print_and_execute(cmd)
    #ImageMath 4 tmp3DImage.nii.gz TimeSeriesDisassemble ops.in_img
    #convert (from ImageMagick)
    #Executables that must exist:
    #ImageMath
    extract_slice = "/home/bkhanal/works/tools/marcoSliceExtractor/myImgSliceExtractor"
    axis = ops.slice_axis
    num = 0
    #print('number of time steps: %s \n' % (str(ops.total_tpts)))
    while num < ops.total_tpts:
        #outputs/results of the executables
        index = str(num+100) #ImageMath extracted slice names start from 100.
        tmp3DImage = '%s%s%s' % (tmp3d_basename, index, file_ext)
        tmp2DImage = 'slice%s.png' % (index,)
        cmd = '%s %s %s %s %s' % (
            extract_slice, tmp3DImage, axis, ops.slice_num, tmp2DImage)
        bu.print_and_execute(cmd, False)
        # Rotate the image for proper orientation.
        if ops.rotate is None:
            cmd = 'convert -rotate 180 %s %s' % (tmp2DImage, tmp2DImage)
        else:
            cmd = 'convert -rotate %s %s %s' % (ops.rotate, tmp2DImage, tmp2DImage)
        bu.print_and_execute(cmd, False)
        if ops.crop:
            cmd = 'convert %s -crop %s +repage %s' % (tmp2DImage, ops.crop, tmp2DImage)
            bu.print_and_execute(cmd, False)
        # Write time-point info
        if ops.time_unit is not None:
            if ops.time_step is not None:
                tpt = float(num) * ops.time_step
            else:
                tpt = num
            cmd = ('convert %s -gravity SouthWest -fill orange -pointsize 12 '
                   '-annotate +0+0 "%s %s" %s' % (
                       tmp2DImage, str(tpt), ops.time_unit, tmp2DImage))
            bu.print_and_execute(cmd, False)

        #Delete individual 3D files.
        bu.print_and_execute('rm ' + tmp3DImage, False)
        #Go to next file
        num += 1

    #Now make the animation and delete individual 2D slices:
    cmd = 'convert slice1*.png -resize %s -set delay %s %s' % (
        ops.resize, ops.delay, ops.out_gif)
    bu.print_and_execute(cmd)
    bu.print_and_execute('rm slice1*.png', False)


if __name__ == "__main__":
    main()
