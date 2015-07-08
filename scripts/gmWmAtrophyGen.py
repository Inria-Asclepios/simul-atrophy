#!/usr/bin/python
""" Generate volume stats, atorphy tables, atrophy maps and/or required
segmentation images from freesurfer results.
This generates results that serve as input to AdLemModel.
"""
# Author: Bishesh Khanal <bishesh.khanal@inria.fr>
#
#
import os
import os.path as op
import subprocess as sp
import argparse as ag

def get_input_options():
    """ Command line interface, get input options and interact with the user.
    """
    parser = ag.ArgumentParser('Generate python data frame '
                               'atrophy maps')
    parser.add_argument('fslstats', help='command for fslstats')
    parser.add_argument('global_mask', help='mask with 1 for tissue and zero '
                        'elsewhere to compute total tissue volume')
    parser.add_argument('label_img', help='label file which has different '
                        'labels of the regions used in atrophy_table file.')
    parser.add_argument('atrophy_in', help='table file that has label ids '
                        'in the first column and the desired contribution '
                        'percentage by this region to the global atrophy')
    parser.add_argument('atrophy_out', help='output table file that will be '
                        'written that has label ids like in the atrophy_in but'
                        ' with actual atrophy values to be prescribed ')
    parser.add_argument('atrophy_img', help='atrophy image that takes regions '
                        'from label_img and values from the created '
                        'atrophy_out table.')
    ops = parser.parse_args()
    return ops


def print_and_execute(command):
    """prints the command and executes it in the shell"""
    print command + '\n'
    sp.call(command, shell=True)


def main():
    """ Get the stats from label image file and compute proper atrophy
    percentage. Run with -h to see options.
    """
    ops = get_input_options()

    a_g = 0 # global atrophy
    alpha_r = {} # Dictionary storing region ids & their contribution to a_g.
    with open(ops.atrophy_in) as fil:
        first_line = (fil.readline()).split()
        if not (first_line[0] == 'global'):
            raise ValueError('First line of atrophy_in file must contain must '
                             'contain a word global followed by the desired '
                             ' global atrophy percentage. ')
        a_g = float(first_line[1])

        for line in fil:
            (region, contrib) = line.split()
            alpha_r[region] = float(contrib)

    # Get the global volume:
    try:
        fsl_out = sp.check_output([ops.fslstats, ops.global_mask, '-v'])
        vol_g = float(fsl_out.split()[0]) # Take the first output of fslstats.
        print('global volume obtained = %f' % (vol_g))
    except sp.CalledProcessError, e:
        print(e.output)

    # Open a file to write atrophy table:
    with open(ops.atrophy_out, 'w') as fil:
        fil.write('labels\tnewValues\n') # the proper first line header.
        for region in alpha_r.keys():
            l_thresh = str(float(region) - 0.5)
            u_thresh = str(float(region) + 0.5)
            try:
                fsl_out = sp.check_output([
                    ops.fslstats, ops.label_img, '-l', l_thresh, '-u',
                    u_thresh, '-v'])
                print('label id %s stats in %s image: ' % (region, ops.label_img))
                print(fsl_out)
                vol_r = float(fsl_out.split()[0])
            except sp.CalledProcessError, e:
                print(e.output)
            a_r = alpha_r[region] * a_g * vol_g / vol_r
            fil.write(region + '\t' + str(a_r) + '\n')

    # Create atrophy map using the created table:
    adlem_dir = os.getenv('ADLEM_DIR')
    if adlem_dir is None:
        raise ValueError('environment variable ADLEM_DIR not set.')
    create_atrophy_map = op.join(adlem_dir,
                                 'build/src/createImageFromLabelImage')
    print_and_execute('%s -t %s -l %s -o %s' %
                      (create_atrophy_map, ops.atrophy_out, ops.label_img,
                       ops.atrophy_img))


if __name__ == "__main__":
    main()
