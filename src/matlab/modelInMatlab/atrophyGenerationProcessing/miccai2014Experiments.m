% Miccai 2014 experiments:
% CSF-Braintissue images images to test our model for different cases.
% Develop masks and atrophy maps.

clear all; clc;
addpath(genpath('/home/bkhanal/Documents/softwares/matlabTools/'));

NBR_LABEL = 0;
CSF_LABEL = 1;
GMWM_LABEL = 2;  
GMWM_ATROPHY = 4;   %3 is for wm if 2 is used only for gm

% CSF sphere in the middle
r_c = 8;

% Brain with atrophy in the middle
r_ba = 4;
r_b = 7;


%External width
w_ext = 3;
% Brain Parenchyma regions: bp:  

% MASKS and their values:
v_c = 1;  v_ba = 4;     v_b = 2;

% The cube size:
xn = 2*r_c + 2*w_ext;
yn = xn;
zn = xn;

% Co-ordinates:
[y x z] = meshgrid(1:yn, 1:xn, 1:zn);       %x is vertically downward!!

% Create CSF masks:
c_c = [w_ext+r_c, w_ext+r_c, w_ext+r_c];

% Before let's be done with atrophy!
atrophy = 0.5;
cBA = (x-c_c(1)).^2 + (y-c_c(2)).^2 + (z-c_c(3)).^2 - r_ba*r_ba;
cBA(cBA<=0) = -1;     cBA(cBA>0) = 0;   cBA(cBA<0) = 1;
cBAimage = cBA .* atrophy;
save('atrophy.mat','cBA');
writemetaimagefile('atrophy.mha',cBAimage,[1 1 1]);

% Get back to CSF mask
cC = (x-c_c(1)).^2 + (y-c_c(2)).^2 + (z-c_c(3)).^2 - r_c*r_c;
cC(cC<=0) = -1;     cC(cC>0) = 0;   cC(cC<0) = 1;

btest3d = cC.*v_c;

% Make everything outside CSF sphere a brain:
btest3d(cC ~= 1) = v_b;

% Create brain masks
cB = (x-c_c(1)).^2 + (y-c_c(2)).^2 + (z-c_c(3)).^2 - r_b*r_b;
cB(cB<=0) = -1;     cB(cB>0) = 0;   cB(cB<0) = 1;
btest3d(cB == 1) = v_b;

save('btest.mat','btest3d');

% Put in the rods:
%% Case 1: A thin rod in the middle
r_r = 1;    w_r = 2*r_c;
c_r = [w_ext+r_c, w_ext+r_c, w_ext+r_c];
cR = ((x-c_r(1))./r_r).^2 + ((z-c_r(3))./r_r).^2 - 1; 
cR(cR<=0) = -1;     cR(cR>0) = 0;   cR(cR<0) = 1;
cR(:,1:w_ext-1,:) = 0;        cR(:,end-w_ext+1:end,:) = 0;

bMask = btest3d;
bMask(cR==1) = v_b;

% Create the image:
imshow3D(bMask)
writemetaimagefile('bMask.mha',bMask,[1 1 1]);

%% Case 2: A thick rod in the middle
r_r = 4;    w_r = 2*r_c;
c_r = [w_ext+r_c, w_ext+r_c, w_ext+r_c];
cR = ((x-c_r(1))./r_r).^2 + ((z-c_r(3))./r_r).^2 - 1; 
cR(cR<=0) = -1;     cR(cR>0) = 0;   cR(cR<0) = 1;
cR(:,1:w_ext-1,:) = 0;        cR(:,end-w_ext+1:end,:) = 0;

bMask = btest3d;
bMask(cR==1) = v_b;

% Create the image:
imshow3D(bMask)
writemetaimagefile('bMask.mha',bMask,[1 1 1]);

%% Case 3: A thick cross rod in the middle
r_r = 5;    w_r = 2*r_c;
c_r = [w_ext+r_c, w_ext+r_c, w_ext+r_c];
cR = ((x-c_r(1))./r_r).^2 + ((z-c_r(3))./r_r).^2 - 1; 
cR(cR<=0) = -1;     cR(cR>0) = 0;   cR(cR<0) = 1;
cR(:,1:w_ext-1,:) = 0;        cR(:,end-w_ext+1:end,:) = 0;

bMask = btest3d;
bMask(cR==1) = v_b;

cR = ((x-c_r(1))./r_r).^2 + ((y-c_r(2))./r_r).^2 - 1; 
cR(cR<=0) = -1;     cR(cR>0) = 0;   cR(cR<0) = 1;
cR(:,:,1:w_ext-1) = 0;        cR(:,:,end-w_ext+1:end) = 0;
bMask(cR==1) = v_b;

cR = ((z-c_r(3))./r_r).^2 + ((y-c_r(2))./r_r).^2 - 1; 
cR(cR<=0) = -1;     cR(cR>0) = 0;   cR(cR<0) = 1;
cR(1:w_ext-1,:,:) = 0;        cR(end-w_ext+1:end,:,:) = 0;
bMask(cR==1) = v_b;

% Create the image:
imshow3D(bMask)
writemetaimagefile('bMask.mha',bMask,[1 1 1]);

%% Case ..
% Central part only brain, csf on two corners of the domain:
csf_w = 5;
bMask(x<=xn) = v_b; %initialize everything with brain value;
bMask((x>2 & x<csf_w+1 & y>2 & y<csf_w+1 ) | (x>2 & x<csf_w+1 & y>yn-csf_w-1 & y<yn-2)) = v_c;
bMask(z<3 | z>zn-3) = v_b;
% bMask(x < 4 | y < 4 | z < 4) = v_c;
imshow3D(bMask);
writemetaimagefile('bMask.mha',bMask,[1 1 1]);

%%  Case .. 
% Central part only brain, csf on four corners of the domain
% everywhere brain except four corners.
csf_w = 5;
bMask(x<=xn) = v_b; %initialize everything with brain value;
bMask((x>2 & x<csf_w+1 & y>2 & y<csf_w+1 ) | (x>2 & x<csf_w+1 & y>yn-csf_w-1 & y<yn-2)...
    | (x>xn-csf_w-1 & x<xn-2 & y>2 & y<csf_w+1) | (x>xn-csf_w-1 & x<xn-2 & y>yn-csf_w-1 & y<yn-2)) = v_c;

bMask(z<3 | z >zn-3) = v_b;
imshow3D(bMask);
writemetaimagefile('bMask.mha',bMask,[1 1 1]);

%% Case ..
% Central part only brain, csf on four corners and a thin layer sorrounding
% the atrophied part of the brain.
csf_w = 5;
bMask(x<=xn) = v_b; %initialize everything with brain value;
bMask((x>2 & x<csf_w+1 & y>2 & y<csf_w+1 ) | (x>2 & x<csf_w+1 & y>yn-csf_w-1 & y<yn-2)...
    | (x>xn-csf_w-1 & x<xn-2 & y>2 & y<csf_w+1) | (x>xn-csf_w-1 & x<xn-2 & y>yn-csf_w-1 & y<yn-2)) = v_c;

bMask(z<3 | z >zn-3) = v_b;

surroundingCSF = (x-c_c(1)).^2 + (y-c_c(2)).^2 + (z-c_c(3)).^2 - (r_ba+1)*(r_ba+1);
surroundingCSF(surroundingCSF<=0) = -1;     surroundingCSF(surroundingCSF>0) = 0;   
surroundingCSF(surroundingCSF<0) = 1;

bMask(surroundingCSF & (~cBA)) = v_c;

% bMask(x < 4 | y < 4 | z < 4) = v_c;
imshow3D(bMask);
writemetaimagefile('bMask.mha',bMask,[1 1 1]);

%% Case ..
% Central part only brain, csf on four corners and a thin layer 
% of CSF on the upper side of the atrophied center
csf_w = 5;
bMask(x<=xn) = v_b; %initialize everything with brain value;
bMask((x>2 & x<csf_w+1 & y>2 & y<csf_w+1 ) | (x>2 & x<csf_w+1 & y>yn-csf_w-1 & y<yn-2)...
    | (x>xn-csf_w-1 & x<xn-2 & y>2 & y<csf_w+1) | (x>xn-csf_w-1 & x<xn-2 & y>yn-csf_w-1 & y<yn-2)) = v_c;

bMask(z<3 | z >zn-3) = v_b;

surroundingCSF = (x-c_c(1)).^2 + (y-c_c(2)).^2 + (z-c_c(3)).^2 - (r_ba+1)*(r_ba+1);
surroundingCSF(surroundingCSF<=0) = -1;     surroundingCSF(surroundingCSF>0) = 0;   
surroundingCSF(surroundingCSF<0) = 1;

bMask(surroundingCSF & (~cBA)) = v_c;
bMask(x>xn/2 | y>yn/2 | z>zn/2) = v_b;
% bMask(x < 4 | y < 4 | z < 4) = v_c;
imshow3D(bMask);
writemetaimagefile('bMask.mha',bMask,[1 1 1]);

%% Case ..
% Central part only brain, csf on four corners and a thin layer 
% of CSF on the lower side of the atrophied center
csf_w = 5;
bMask(x<=xn) = v_b; %initialize everything with brain value;
bMask((x>2 & x<csf_w+1 & y>2 & y<csf_w+1 ) | (x>2 & x<csf_w+1 & y>yn-csf_w-1 & y<yn-2)...
    | (x>xn-csf_w-1 & x<xn-2 & y>2 & y<csf_w+1) | (x>xn-csf_w-1 & x<xn-2 & y>yn-csf_w-1 & y<yn-2)) = v_c;

bMask(z<3 | z >zn-3) = v_b;

surroundingCSF = (x-c_c(1)).^2 + (y-c_c(2)).^2 + (z-c_c(3)).^2 - (r_ba+1)*(r_ba+1);
surroundingCSF(surroundingCSF<=0) = -1;     surroundingCSF(surroundingCSF>0) = 0;   
surroundingCSF(surroundingCSF<0) = 1;

bMask(surroundingCSF & (~cBA) & x<xn/2 & y<yn/2 & z<zn/2) = v_c;
imshow3D(bMask);
writemetaimagefile('bMask.mha',bMask,[1 1 1]);
%% Generate image to be warped to see the effect of the velocity field:
img = zeros(xn,yn,zn);
two = 3*ones(xn,yn,zn);
img(logical(rem(x,two))) = 1;
img(logical(rem(y,two))) = 1;
% img(logical(rem(z,two))) = 1;
imshow3D(img);
writemetaimagefile('checkerImage.mha',img,[1 1 1]);
