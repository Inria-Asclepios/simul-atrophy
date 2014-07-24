% Extracts and saves a .mat file with a 2D image slice, corresponding
% segmentation mask and divergence map from following inputs:
% 
% 1. path (which must have three files of following names:
% baseline.nii, baseline_seg.nii and divergence.nii) this is where a .mat
% file is saved with a segAndDiv prefix and the info of the following two
% inputs.
% 
% 2. cutting plane that is 'coronal' or 'sagittal' or 'axial'
% 
% 3. slice number that you want to extract
% 
% 4. (optional) DISPLAY if passed as true, all the extracted will be shown
% 


function [img_slice seg_mask div_slice] = getSegAndDivMask(dirname,cutting_plane,slice_num,DISPLAY)
if (nargin < 4)
    DISPLAY = false;
end
% dirname = 'RealData/Template';

img_nii = load_nii([dirname '/baseline.nii']);
seg_nii = load_nii([dirname '/baseline_seg.nii']);
div_nii = load_nii([dirname '/divergence.nii']);

% cutting_plane = 'coronal';
% slice_num = 110;

[img_img3D img_slice] = getImageSliceFromNifti(img_nii,cutting_plane,slice_num);
[seg_img3D seg_slice] = getImageSliceFromNifti(seg_nii,cutting_plane,slice_num);
[div_img3D div_slice] = getImageSliceFromNifti(div_nii,cutting_plane,slice_num);

seg_mask = false(size(seg_slice));
seg_mask((seg_slice == 2) | (seg_slice == 3)) = true;

if (DISPLAY == true)
    subplot(221), imshow(img_slice,[]), title('baseline image slice');
    subplot(222), imshow(seg_slice,[]), title('baseline image segmentation');
    subplot(223), imshow(div_slice,[]), title('divergence map from registration');
    subplot(224), imshow(seg_mask,[]), title('mask for CSF-parenchyma');
end


end