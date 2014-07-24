% MUST ADAPT FOR THE ORIENTATION!!! FROM THE HEADER GET ORIENTATION INFO
% AND CHANGE TO THE STANDARD ONE!
% FOR STORING THREE PLANES!!! OR RELATE THIS WITH THE MATLAB THINGY!!
% ------------------------------------------------------------------------
% Get a 3D image with rearranged axes [x y z] for the chose plane from a
% given 3D MRI image in Nifti format; Optionally get a 2D slice image too.
% *************************
% Input descriptions:
% -----------------------
% mri_nii : a Nifti format data. e.g using one obtained by using function
% load_nii
% -----------------------
% plane: 'axial' or 'coronal' or 'sagittal'.
% -----------------------
% slice_num (optional): The slice number to be extracted.
% ************************
% Output descriptions:
% img: 3D image
% 
% Bishesh Khanal
% Asclepios, INRIA
% ----------------------

function [img img_slice] = getImageSliceFromNifti(mri_nii,plane,slice_num)

% mr_nii = load_nii('../Tissues_Images/MRI.nii');
mr_img = mri_nii.img;
switch lower(plane) %change to lower case
    case 'coronal'
        mr_coronal = permute(mr_img,[3 1 2]);
        mr_coronal = flipdim(mr_coronal,1);
        mr_coronal = flipdim(mr_coronal,2);
        img= mr_coronal;

    case 'axial'
        mr_axial = permute(mr_img,[2 1 3]);
        mr_axial = flipdim(mr_axial,1);
        mr_axial = flipdim(mr_axial,2);
        img = mr_axial;

    case 'sagittal'
        mr_sagittal = permute(mr_img,[3 2 1]);
        mr_sagittal = flipdim(mr_sagittal,1);
        mr_sagittal = flipdim(mr_sagittal,2);
        img = mr_sagittal;

    otherwise
        error('invalide plane option');
end

if (nargin == 3)
    img_slice = img(:,:,slice_num);
end

end

