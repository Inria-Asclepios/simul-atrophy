#!/usr/bin/python
""" Generate different kinds of synthetic images.

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
    parser = ag.ArgumentParser('''Create synthetic images to create different
    atrophy patterns.
    ''')
    parser.add_argument('img_size', help='img_size in the format mxnxr')
    parser.add_argument('out_file', help='Output image file')
    ops = parser.parse_args()
    return ops


def box_with_spheres(r1, r2, r3):
    """ Return a numpy 3D array with three concentric spheres in it.
    r1 < r2 < r3.
    The dimension of the box is a cube with length = d3 + 1.
    """
    xx, yy, zz = np.mgrid[:2*r3, :2*r3, :2*r3]

    # A sphere
    cc = r3
    sph = (xx-cc)**2 + (yy-cc)**2 + (zz-cc)**2

    sph1 = sph <= r1**2
    ring1 = (sph > r1**2) & (sph <= r2**2)
    ring2 = (sph > r2**2) & (sph <= r3**2)

    sph = (50 * sph1.astype(float) + 100 * ring1.astype(float) +
           10 * ring2.astype(float))
    # sph = (2 * sph1.astype(float) + 3 * ring1.astype(float) +
    #        1 * ring2.astype(float))

    return sph

def box_with_tube(r1, r2, r3):
    """ Return a numpy 3D array with three concentric cylinders in it."""
    h = r3*2 + 2
    xx, yy, zz = np.mgrid[:h, :h, :h] # axes: 0 to h-1 

    # A cylinder
    cc = r3
    cyl_fn = (xx-cc)**2 + (yy-cc)**2
    #h1 = (zz > 4) & (zz < (h-5))
    h1 = (zz > 6) & (zz < (h-7))
    h2 = ((zz <= 6 ) & (zz >= 3)) | ((zz >= h-7) & (zz <= h-4))
    #h2 = (zz ==  4) | (zz == 3) | (zz == h-5) | (zz == h-4)
    h3 = (zz ==  2) | (zz == 1) | (zz == h-3) | (zz == h-2)

    cyl1 = (cyl_fn <= r1**2) & h1
    cyl2 = (cyl_fn > r1**2) & (cyl_fn <= r2**2 ) & h1
    cyl3 = (cyl_fn > r2**2) & (cyl_fn <= r3**2) & (h1 | h2)

    extra = (cyl_fn <= r2**2) & h2

    cyl = (50 * cyl1.astype(float) + 100 * cyl2.astype(float) +
           10 * cyl3.astype(float) + 10 * h3.astype(float) +
           70 * extra.astype(float))

    # cyl = (2 * cyl1.astype(float) + 3 * cyl2.astype(float) +
    #        1 * cyl3.astype(float) + 1 * h3.astype(float) +
    #        2 * extra.astype(float))

    return cyl



def main():
    """ Launch recon_all long jobs. Run with -h to see options."""
    ops = get_input_options()
    #img_data = box_with_spheres(5, 10, 14)
    img_data = box_with_tube(5, 10, 14)

    # Create and save the image
    img = nib.Nifti1Image(img_data, np.eye(4))
    nib.save(img, ops.out_file)





if __name__ == "__main__":
    main()
