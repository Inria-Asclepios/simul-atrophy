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
writemetaimagefile('cyl.mha',cyl,[1 1 1]);

%% Some piecewise uniform atrophy maps:
% Uniform in individual layers but different accross different layers.
atrophy1 = cyl/10;
writemetaimagefile('aL.mha',atrophy1,[1 1 1]);
% Uniform across all layers
atrophy2 = 0.1*ones(size(atrophy1));
writemetaimagefile('aU.mha',atrophy2,[1 1 1]);
%% Gaussian atrophy maps:
sigma = 4;
max_atrophy = 0.4;
center = [0 0 0];
atrophy1 = max_atrophy*exp(-((x-center(1)).^2 + (y-center(2)).^2 + ...
    (z-center(3)).^2)/(2*sigma*sigma));
writemetaimagefile('aGX0Y0Z0S4.mha',atrophy1,[1 1 1]);
%% Now let's convert one of the layers into CSF.
% Case 1: Outermost layer as CSF:
mask5asCsf = cyl;
mask5asCsf(cyl==5) = 1;
mask5asCsf(cyl==3 | cyl==4) = 2;
writemetaimagefile('m5.mha',mask5asCsf,[1 1 1]);

% Case 2: Innermost layer as CSF:
mask2asCsf = cyl;
mask2asCsf(cyl==2) = 1;
mask2asCsf(cyl==3 | cyl==4 | cyl==5) = 2;
writemetaimagefile('m2.mha',mask2asCsf,[1 1 1]);

% Case 3: Second innermost layer as CSF:
mask3asCsf = cyl;
mask3asCsf(cyl==3) = 1;
mask3asCsf(cyl==2 | cyl==4 | cyl==5) = 2;
writemetaimagefile('m3.mha',mask3asCsf,[1 1 1]);

%% Let's have isolated CSF spheres inside the cylinder.
r = (zn/4)-1;
c1 = [0 0 zn/4];
c2 = [0 0 -zn/4];
disjointSpheres = zeros(xn,yn,zn);
disjointSpheres = (((x-c1(1)).^2 + (y-c1(2)).^2 + (z-c1(3)).^2) < r^2) | ...
    (((x-c2(1)).^2 + (y-c2(2)).^2 + (z-c2(3)).^2) < r^2) + disjointSpheres; 
% imshow3D(disjointSpheres);
writemetaimagefile('twoSpheres.mha',double(disjointSpheres),[1 1 1]);
mask2SpheresAsCsf = cyl;
mask2SpheresAsCsf(disjointSpheres==1) = 1;
writemetaimagefile('mS2.mha',double(mask2SpheresAsCsf),[1 1 1]);
%% Eigen vectors and then tensor images.
u_field = zeros(xn,yn,zn,3);

% Principle direction along the length of the cylinder: (pointing z-axis);
u_field(:,:,:,1) = 0;   % for radial inwards -y;
u_field(:,:,:,2) = 1;   % for radial inwards -x;
u_field(:,:,:,3) = 0;   % for radial inwards 0;
u_field = permute(u_field,[4 1 2 3]); % to save it as vector field.
%  Rotating vector fields as desired. 
% u_rotated = rotate3dVectorField(u_field,pi/2,pi/2,pi,[1 2 3]);
% writemetaimagefile('vec001test.mha',u_rotated,[1 1 1],[0 0 0],3);

% Get two orthogonal unit vector from the given input vector.
u_unit = getUnitVectors(u_field);
% writemetaimagefile('vecUpwards.mha',u_unit,[1 1 1],[0 0 0],3);

%get two orthogonal unit vectors to the given input vector:
[v w] = getTwoOrthogonalVectorFields(u_unit);
% writemetaimagefile('vec001v.mha',v,[1 1 1],[0 0 0],3);
% writemetaimagefile('vec001w.mha',w,[1 1 1],[0 0 0],3);

% create tensor images from three eigenvectors and three eigenvalues:
tensorImage = getDiffusionImage(u_unit,v,w,10,4,1);
tensorImageForNii = zeros([size(tensorImage) 1]);
% clear tensorImageForNii;
tensorImageForNii(:,:,:,:,1) = tensorImage;
tensorImageForNii = permute(tensorImageForNii,[2 3 4 5 1]);
tensor_nii = make_nii(tensorImageForNii);
tensor_nii.img = tensorImageForNii;
tensor_nii.hdr.dime.dim = [5 size(tensorImageForNii) 1 1];
tensor_nii.hdr.dime.intent_code = 1007;
save_nii(tensor_nii,'tFxL10L4L1.nii');
% save_nii(tensor_nii,'tRiL10L4L1.nii');
% save_nii(tensor_nii,'tUzL10L4L1.nii');

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


    