#!/usr/bin/python
""" Compare Jacobian determinant of registartion algorithm with the prescribed
one.

"""
# Author: Bishesh Khanal <bishesh.khanal@inria.fr>
#
#

import argparse as ag
import os
import os.path as op
import numpy as np
import nibabel as nib


def get_input_options():
    """ Command line interface, get input options and interact with the user.
    """
    parser = ag.ArgumentParser('''Compare registration vs prescribed simulated
    atrophy.
    ''')
    parser.add_argument('fixed_img', help='fixed image, the baseline.')
    parser.add_argument('moving_img', help='moving image, simulated')
    ops = parser.parse_args()
    return ops


def main():
    """ Launch recon_all long jobs. Run with -h to see options."""
    ops = get_input_options()
    #img_data = box_with_spheres(5, 10, 14)





if __name__ == "__main__":
    main()
