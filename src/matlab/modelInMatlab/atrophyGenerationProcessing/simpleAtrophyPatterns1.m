% simpleAtrophyPatterns1.m
% Create a certain pattern: (See towards the end for 2D)

clear all; clc;
addpath(genpath('/home/bkhanal/Documents/softwares/matlabTools/'));

% borderWidth:
bw = 2;

% Left: L Right: R Middle: M
% CSF: c:
% Two spheres of CSF on left and right
r_cL = 8;      r_cR = 8;       
% A cylinder joining the two above spheres
r_cM = 1;    w_cM = 4;

% Brain Parenchyma regions: bp:  
% Two spheres concentric to the CSF spheres but of smaller size:
r_bpL = 5;     r_bpR = 5;

% Another Brain Parenchyma region: bp1
r_bp1L = 2;     r_bp1R = 2;

% MASKS and their values:
% v_cL = 20;  v_cM = 40;  v_cR = 60;  v_bpL = 100;    v_bpR = 100;
v_cL = 1;  v_cM = 1;  v_cR = 1;  v_bpL = 2;    v_bpR = 2;
% v_bp1L = 200;           v_bp1R = 200;
v_bp1L = 3;           v_bp1R = 3;

% The cube size:
xn = 2*r_cL + 2*bw;
yn = 2*r_cL + w_cM + 2*r_cR + 2*bw;
zn = xn;

% Co-ordinates:
[y x z] = meshgrid(1:yn, 1:xn, 1:zn);       %x is vertically downward!!


% Create CSF masks:
c_cL = [bw+r_cL, bw+r_cL, bw+r_cL];     c_cR = [bw+r_cR, yn-(bw+r_cR), bw+r_cR];
c_cM = [bw+r_cL, bw + r_cL*2 + w_cM/2, bw+r_cL];

cL = (x-c_cL(1)).^2 + (y-c_cL(2)).^2 + (z-c_cL(3)).^2 - r_cL*r_cL;
cL(cL<=0) = -1;     cL(cL>0) = 0;   cL(cL<0) = 1;


cR = (x-c_cR(1)).^2 + (y-c_cR(2)).^2 + (z-c_cR(3)).^2 - r_cR*r_cR;
cR(cR<=0) = -1;     cR(cR>0) = 0;   cR(cR<0) = 1;

cM = ((x-c_cM(1))./r_cM).^2 + ((z-c_cM(3))./r_cM).^2 - 1; 
cM(cM<=0) = -1;     cM(cM>0) = 0;   cM(cM<0) = 1;
cM(:,1:bw+2*r_cL-1,:) = 0;        cM(:,end-bw-2*r_cR+1:end,:) = 0;

btest3d = v_cL.*cL + v_cR.*cR;
btest3d(cM==1) = v_cM;

% Create bp masks
c_bpL = c_cL;       c_bpR = c_cR;
bpL = (x-c_bpL(1)).^2 + (y-c_bpL(2)).^2 + (z-c_bpL(3)).^2 - r_bpL*r_bpL;
bpL(bpL<=0) = -1;     bpL(bpL>0) = 0;   bpL(bpL<0) = 1;

bpR = (x-c_bpR(1)).^2 + (y-c_bpR(2)).^2 + (z-c_bpR(3)).^2 - r_bpR*r_bpR;
bpR(bpR<=0) = -1;     bpR(bpR>0) = 0;   bpR(bpR<0) = 1;

btest3d(bpL==1) = v_bpL;     btest3d(bpR==1) = v_bpR;

c_bp1L = c_cL;       c_bp1R = c_cR;
bp1L = (x-c_bp1L(1)).^2 + (y-c_bp1L(2)).^2 + (z-c_bp1L(3)).^2 - r_bp1L*r_bp1L;
bp1L(bp1L<=0) = -1;     bp1L(bp1L>0) = 0;   bp1L(bp1L<0) = 1;

bp1R = (x-c_bp1R(1)).^2 + (y-c_bp1R(2)).^2 + (z-c_bp1R(3)).^2 - r_bp1R*r_bp1R;
bp1R(bp1R<=0) = -1;     bp1R(bp1R>0) = 0;   bp1R(bp1R<0) = 1;

btest3d(bp1L==1) = v_bp1L;     btest3d(bp1R==1) = v_bp1R;

% Create the image:

imshow3D(btest3d)
writemetaimagefile('btest.mha',btest3d,[1 1 1]);
% btest3dnii = make_nii(btest3d);
save_nii(make_nii(btest3d),'btest.nii');
%% Save the image:
save('btest.mat','btest3d');
writemetaimagefile('btest.mha',btest3d,[1 1 1]);

b = btest3d;
% atrophy:
% Put atrophy everywhere in the brain, and compute it's sum.
a_bpL = 0.2;        a_bpR = 0.2;
a = zeros(xn,yn,zn);
a(b==v_bpL) = a_bpL;
a(b==v_bpR) = a_bpR;
aInBrain = sum(a(:));

% CASE 1:
% Uniformly distribute atrophy in the csf:
a1 = zeros(xn,yn,zn);
a1(b==v_cL | b==v_cM | b==v_cR) = 1;
a1(a1==1) = -aInBrain/sum(a1(:));
a1(b==v_bpL) = a_bpL;
a1(b==v_bpR) = a_bpR;

% Check that total atrophy should be zero:
display('a1 sum');
sum(a1(:))
writemetaimagefile('testAtrophy1.mha',a1,[1 1 1]);

% CASE 2:
% Zero Exapnsion on left:
a2 = zeros(xn,yn,zn);
a2(b==v_cM | b==v_cR) = 1;
a2(a2==1) = -aInBrain/sum(a2(:));
a2(b==v_bpL) = a_bpL;
a2(b==v_bpR) = a_bpR;
display('a2 sum:');
sum(a2(:))
writemetaimagefile('testAtrophy2.mha',a2,[1 1 1]);

% CASE 3:
% Zero Expansion on right:
a3 = zeros(xn,yn,zn);
a3(b==v_cL | b==v_cR) = 1;
a3(a3==1) = -aInBrain/sum(a3(:));
a3(b==v_bpL) = a_bpL;
a3(b==v_bpR) = a_bpR;
display('a3 sum:');
sum(a3(:))
writemetaimagefile('testAtrophy3.mha',a3,[1 1 1]);

% CASE 4:
% Zero Expansion on both sides:
a4 = zeros(xn,yn,zn);
a4(b==v_cM) = 1;
a4(a4==1) = -aInBrain/sum(a4(:));
a4(b==v_bpL) = a_bpL;
a4(b==v_bpR) = a_bpR;
display('a4 sum:');
sum(a4(:))
writemetaimagefile('testAtrophy4.mha',a4,[1 1 1]);

%%  ******************** In 2D: ***********************************

slice = zn/2;
btest2d = btest3d(:,:,slice); btest2d = padarray(btest2d,[1 1],'pre');
bpL = bpL(:,:,slice);   bpL = padarray(bpL,[1 1],'pre');
bpR = bpR(:,:,slice);   bpR = padarray(bpR,[1 1],'pre');
cL = cL(:,:,slice);     cL = padarray(cL,[1 1],'pre');
cM = cM(:,:,slice);     cM = padarray(cM,[1 1],'pre');
cR = cR(:,:,slice);     cR = padarray(cR,[1 1],'pre');
xn = xn+1; yn = yn+1;
%% Put Atrophies
a = zeros(xn,yn); 
b = btest2d;
bMask = zeros(size(btest2d));     bMask(b==v_cL | b==v_cM | b==v_cR) = 1;
bMask(b==v_bpL | b==v_bpR) = 2;

% Put atrophy everywhere in the brain, and compute it's sum.
a_bpL = 0.2;        a_bpR = 0.2;
a(b==v_bpL) = a_bpL;
a(b==v_bpR) = a_bpR;
aInBrain = sum(a(:));

% CASE 1:
% Uniformly distribute atrophy in the csf:
a1 = zeros(xn,yn);
a1(b==v_cL | b==v_cM | b==v_cR) = 1;
a1(a1==1) = -aInBrain/sum(a1(:));
a1(b==v_bpL) = a_bpL;
a1(b==v_bpR) = a_bpR;

% Check that total atrophy should be zero:
display('a1 sum');
sum(a1(:))

% CASE 2:
% Zero Exapnsion on left:
a2 = zeros(xn,yn);
a2(b==v_cM | b==v_cR) = 1;
a2(a2==1) = -aInBrain/sum(a2(:));
a2(b==v_bpL) = a_bpL;
a2(b==v_bpR) = a_bpR;
display('a2 sum:');
sum(a2(:))

% CASE 3:
% Zero Expansion on right:
a3 = zeros(xn,yn);
a3(b==v_cL | b==v_cR) = 1;
a3(a3==1) = -aInBrain/sum(a3(:));
a3(b==v_bpL) = a_bpL;
a3(b==v_bpR) = a_bpR;
display('a3 sum:');
sum(a3(:))

% CASE 4:
% Zero Expansion on both sides:
a4 = zeros(xn,yn);
a4(b==v_cM) = 1;
a4(a4==1) = -aInBrain/sum(a4(:));
a4(b==v_bpL) = a_bpL;
a4(b==v_bpR) = a_bpR;
display('a4 sum:');
sum(a4(:))

%% Let's see:
minC = min([min(a1(:)),min(a2(:)),min(a3(:)),min(a4(:))]);
maxC = max([max(a1(:)),max(a2(:)),max(a3(:)),max(a4(:))]);
clim = [minC maxC];
figure;
subplot(221), imagesc(a1,clim), axis image;
subplot(222), imagesc(a2,clim), axis image;
subplot(223), imagesc(a3,clim), axis image;
subplot(224), imagesc(a4,clim), axis image;

%% Let's save:
save('btest2d.mat','btest2d','bMask','a1','a2','a3','a4');