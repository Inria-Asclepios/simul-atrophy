'''
image_math.py
Some helper functions to do some maths with 3D images.

Author:
Bishesh Khanal
Asclepios, Inria Sophia Antipolis
'''

import numpy as np
import nibabel as ni
import subprocess as sp
import shlex # To convert command string to list which is required when using
# shell = False SEE: http://pymotw.com/2/shlex/ for more of it.
import os
import os.path as op
import time as tm

def get_min_max(imgs, are_fnames=True):
    '''
    Return the minimum and maximum intensities from among the given images.
    imgs must be a list!
    are_fnames: boolean; if True imgs is a list of img_filenames.
    So first loads those files, otherwise considers them as numpy array images.
    '''
    mins, maxs = np.zeros(len(imgs)), np.zeros(len(imgs))
    i = 0
    for img in imgs:
        if are_fnames:
            img = (ni.load(img)).get_data()
        mins[i], maxs[i] = np.amin(img), np.amax(img)
        i += 1
    return np.amin(mins), np.amax(maxs)


def diff_image(in_img1, in_img2, out_img=None):
    '''
    If third argument is not given or None (default), considers
    first two inputs as numpy nd array, so I simply
    return in_img1 - in_img2.

    Else, if third argument given, assumes that all input arguments are
    filenames.
    Write out_img to disk.
    out_image = in_img1 - in_img2
    Copies header info from in_img1.
    '''
    if out_img is None:
        img3_data = in_img1 - in_img2
        return img3_data
    else:
        img1 = ni.load(in_img1)
        img1_data = img1.get_data()

        img2 = ni.load(in_img2)
        img2_data = img2.get_data()

        img3_data = img1_data - img2_data
        img3 = ni.Nifti1Image(img3_data, img1.get_affine(), img1.get_header())
        ni.save(img3, out_img)
        return

def mult_scalar_imgs(im_fnames, out_im_fnames, mult_by):
    '''
    Multiply input images by the given scalar and save them.
    '''
    for imf, imf_out in zip(im_fnames, out_im_fnames):
        img = ni.load(imf)
        out_data = img.get_data() * mult_by
        # Create new output image
        img_out = ni.Nifti1Image(out_data, img.get_affine(), img.get_header())
        ni.save(img_out, imf_out)


def reciprocal_imgs(im_fnames, rec_im_fnames):
    '''
    Take reciprocal of every values in the images and save the resulting images.
    '''
    for imf, imf_rec in zip(im_fnames, rec_im_fnames):
        img = ni.load(imf)
        rec_data = np.reciprocal(img.get_data())
        # Create new output image image
        img_rec = ni.Nifti1Image(rec_data, img.get_affine(), img.get_header())
        ni.save(img_rec, imf_rec)

def create_jac_dets(vel_fnames, jac_fnames, ants_dir, log_jac=False):
    '''
    Create Log/Jacobian determinants of the given vector_fields and save them to
    the given jac_fnames.
    <vel_fnames> [] of vectorfield filenames
    <jac_fnames> [] of desired filenames for results.
    <ants_dir> ANTS bin directory that contains the executables.
    <log_jac> bool to decide log Jacobian determinants or normal one.
    '''
    create_jac = op.join(ants_dir, 'CreateJacobianDeterminantImage')
    for vel, jac in zip(vel_fnames, jac_fnames):
        #print op.basename(vel), op.basename(jac)
        cmd = '%s 3 %s %s %i' % ( #%i since ants command accepts 0 or 1.
            create_jac, vel, jac, log_jac)
        print cmd
        sp.call(shlex.split(cmd), shell=False) #shell=False needs list
        tm.sleep(0.1)

