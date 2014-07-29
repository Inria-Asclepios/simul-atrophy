% First create layered cylindrical volumes with different labels.
% Implicit form of a cylinder:
% x^2 + y^2 <= r^2, -h/2< <= z <= h/2

% clear all; clc;
addpath(genpath('/home/bkhanal/Documents/softwares/matlabTools/'));

% Have different layers within a cylinder each with a different label.
% Total number of layers
numOfLayers = 4;
labels = 2:numOfLayers+1;  % label 0 is for non-brain and 1 for csf.
% width of each layer
layerWidth = 3;

% maximum radius is therefore:
max_r = numOfLayers*layerWidth;
% image size, even dimensions:
xn = max_r*2 + 4; 
yn = xn; zn = xn;
% height of a cylinder:
h = xn-4;

[x y z] = meshgrid(-xn/2:xn/2, -yn/2:yn/2, -zn/2:zn/2);
xn = xn+1; yn = yn+1; zn = zn+1;

cyl = zeros(xn,yn,zn);
cylOld = cyl;
for index = 1:numOfLayers %:-1:1
    r = index*layerWidth;   % radius of the cylinder.
    cylTmp  = (x.*x + y.*y <= r*r) & (z > -h/2 & z < h/2);
    cyl = (labels(index) * (cylTmp & ~cylOld)) + cyl;
    cylOld = cylTmp;
end

% Let's start with a single cylinder!!
% imshow3D(cyl)
writemetaimagefile('cylinder.mha',cyl,[1 1 1]);

%% Some atrophy maps:
% Uniform in individual layers but different accross different layers.
atrophy1 = cyl/10;
writemetaimagefile('atrophy1.mha',atrophy1,[1 1 1]);

%% Now let's convert one of the layers into CSF.
% Case 1: Outermost layer as CSF:
mask5asCsf = cyl;
mask5asCsf(cyl==5) = 1;
mask5asCsf(cyl==3 | cyl==4) = 2;
writemetaimagefile('mask5asCsf.mha',mask5asCsf,[1 1 1]);

% Case 2: Innermost layer as CSF:
mask2asCsf = cyl;
mask2asCsf(cyl==2) = 1;
mask2asCsf(cyl==3 | cyl==4 | cyl==5) = 2;
writemetaimagefile('mask2asCsf.mha',mask2asCsf,[1 1 1]);

% Case 3: Second innermost layer as CSF:
mask3asCsf = cyl;
mask3asCsf(cyl==3) = 1;
mask3asCsf(cyl==2 | cyl==4 | cyl==5) = 2;
writemetaimagefile('mask3asCsf.mha',mask3asCsf,[1 1 1]);

%% Let's have isolated CSF spheres inside the cylinder.

%% Eigen vectors
u_field = zeros(xn,yn,zn,3);

% Principle direction along the length of the cylinder: (pointing z-axis);
u_field(:,:,:,1) = 0;   % for radial inwards -y;
u_field(:,:,:,2) = 0;   % for radial inwards -x;
u_field(:,:,:,3) = 1;
u_field = permute(u_field,[4 1 2 3]); % to save it as vector field.
%  Rotating vector fields as desired. 
% u_rotated = rotate3dVectorField(u_field,pi/2,pi/2,pi,[1 2 3]);
% writemetaimagefile('vec001test.mha',u_rotated,[1 1 1],[0 0 0],3);

% Get two orthogonal unit vector from the given input vector.
u_unit = getUnitVectors(u_field);
writemetaimagefile('vecUpwards.mha',u_unit,[1 1 1],[0 0 0],3);

%get two orthogonal unit vectors to the given input vector:
[v w] = getTwoOrthogonalVectorFields(u_unit);
writemetaimagefile('vec001v.mha',v,[1 1 1],[0 0 0],3);
writemetaimagefile('vec001w.mha',w,[1 1 1],[0 0 0],3);

%% create tensor images from three eigenvectors and three eigenvalues:
tensorImage = getDiffusionImage(u_unit,v,w,10,7,5);
tensorImageForNii = zeros([size(tensorImage) 1]);
% clear tensorImageForNii;
tensorImageForNii(:,:,:,:,1) = tensorImage;
tensorImageForNii = permute(tensorImageForNii,[2 3 4 5 1]);
tensor_nii = make_nii(tensorImageForNii);
tensor_nii.img = tensorImageForNii;
tensor_nii.hdr.dime.dim = [5 size(tensorImageForNii) 1 1];
tensor_nii.hdr.dime.intent_code = 1007;
% save_nii(tensor_nii,'tensorRadialInwards.nii');
save_nii(tensor_nii,'tensorUpwardsZ.nii');

%% Vector field 
% writemetaimagefile('tensorImage.mha',tensorImage,[1 1 1],[0 0 0],6);
u = zeros(xn,yn,zn,1,3);
% Principle direction along the length of the cylinder: (pointing z-axis);
u(:,:,:,:,3) = 10;
u_nii = make_nii(u);
u_nii.img = u;
u_nii.hdr.dime.dim = [5 size(u) 1 1];
u_nii.hdr.dime.intent_code = 1007;
save_nii(u_nii,'vec001test.nii');


    