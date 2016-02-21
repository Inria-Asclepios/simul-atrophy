#!/usr/bin/env python
import os
import os.path as op
import argparse as ag
import bish_utils as bu
import time as tm

def get_input_options():
    """ command line interface, get input options and interact with the user
    """
    #parse the command line inputs
    parser = ag.ArgumentParser(
        'Read list of input big_images and patient_ids from a file and paste'
        'cropped simulated resulting images to the input original image ')
    parser.add_argument(
        'imgs_simIds', help='File with two columns. 1st: images of bigger size'
        ' with full path. 2nd: sim_dir where the results of simulated images '
        ' of cropped smaller size are present in $ADLEM_DIR/results/patients/')
    parser.add_argument('crop_mask', help='File name of the mask that was used '
                        'to crop original images. ')
    parser.add_argument(
        'sim_img', help='File name of the simluated image of smaller size.')
    parser.add_argument(
        'pasted_img', help='Output file name of the pasted images.')
    parser.add_argument(
        '--overwrite', action='store_true', help='overwrites output if already '
        'present. Otherwise skip pasting if pasted_img exists.')
    parser.add_argument('--print_only', action='store_true', help='just print'
                        ', dont do anything!')
    ops = parser.parse_args()
    return ops


def main():
    """ Runs the script. Run with -h to see options.
    """
    ops = get_input_options()
    adlem_dir = os.getenv('ADLEM_DIR')
    if adlem_dir is None:
        raise ValueError('environment variable ADLEM_DIR not set.')
    res_dir = op.join(adlem_dir, 'results/patients')
    paste_img = op.join(adlem_dir, 'build/src/pasteImageToBiggerImage')
    lines = bu.get_lines_as_list(ops.imgs_simIds, '#')
    for line in lines:
        rb, sim_id = line.split()[0], line.split()[1]
        cropped_sf = op.join(res_dir, sim_id, ops.sim_img)
        pasted_sf = op.join(res_dir, sim_id, ops.pasted_img)
        crop_mask = op.join(res_dir, sim_id, ops.crop_mask)
        if op.exists(cropped_sf) and op.exists(rb) and op.exists(crop_mask):
            if not op.exists(pasted_sf) or ops.overwrite:
                cmd = ('%s %s %s %s %s\n\n'
                       % (paste_img, cropped_sf, rb, crop_mask, pasted_sf))
                if ops.print_only:
                    print cmd
                else:
                    bu.print_and_execute(cmd)
                    tm.sleep(0.1)

if __name__ == "__main__":
    main()
