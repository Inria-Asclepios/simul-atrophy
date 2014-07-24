
% Load Image data:
clear all;
addpath(genpath('/home/bkhanal/works/AdLemModel/src/matlab'));
dirname = '/home/bkhanal/works/AdLemModel/data/Template';
view_dir = 'AP'; %viewing direction: anterior-posterior


img_nii = load_nii([dirname '/baseline.nii']);
seg_nii = load_nii([dirname '/baseline_seg.nii']);
div_nii = load_nii([dirname '/divergence.nii']);

img = getImageSliceFromNifti(img_nii,'coronal');
seg_img = getImageSliceFromNifti(seg_nii,'coronal');
div_map = getImageSliceFromNifti(div_nii,'coronal');

seg_mask = false(size(seg_img));
seg_mask((seg_img == 2) | (seg_img == 3)) = true;
% load(savedFile);
USE_REAL_DATA = 0;


[brain_mask min_max_yxz] = getCuboidConvexHull(seg_mask);
[yn xn zn] = size(brain_mask);
img_cr = double(img(min_max_yxz(1) : min_max_yxz(2),...
    min_max_yxz(3) : min_max_yxz(4), min_max_yxz(5) : min_max_yxz(6)));

%pad with extra_CSF width of zeros all around.
eCSF = 2;
img_cr = padarray(img_cr,[eCSF eCSF eCSF]);

 
    %Grid nodes have one dimension greater than the cell centers!
    yn = yn+eCSF*2+1;
    xn = xn+eCSF*2+1;
    zn = zn+eCSF*2+1;
    
% add a dummy
img_cr = cat(3,zeros(yn-1,xn-1,1),img_cr);
img_cr = [zeros(1,xn,zn); zeros(yn-1,1,zn) img_cr];


muBrain = 2500;
lambdaBrain = 2000;

% muB/muc and lambdaB/lambdaC
muRatio = 1; lambdaRatio = 1;
% muRatio = 1000; lambdaRatio = 1000;
% muRatio = 2; lambdaRatio = 2;

% staggeredAugLagr3D(muBrain,muRatio,lambdaBrain,lambdaRatio,seg_mask,div_map,USE_REAL_DATA);

% AdLemTaras3D(seg_mask,div_map,USE_REAL_DATA,muBrain,muRatio,lambdaBrain,lambdaRatio,voi);

% [v2 p] = staggeredVar2Dfaster(seg_mask,div_slice_m,USE_REAL_DATA,muBrain,muRatio,lambdaBrain,lambdaRatio);


% imshow(seg_mask, []);

% saved_file = [dirname '/' cutting_plane num2str(slice_num) '_muR' num2str(muRatio) '.mat'];
% save(saved_file, 'img_slice', 'seg_mask','div_slice', 'v2', 'p');
% 
% % Load v2 and img_slice of certain muRatio:
% clear v2; clear p;
% muRatio = 1000; 
% saved_file = [dirname '/' cutting_plane num2str(slice_num) '_muR' num2str(muRatio) '.mat'];
% load(saved_file);
% 
