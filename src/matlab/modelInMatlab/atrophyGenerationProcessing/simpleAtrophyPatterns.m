% simpleAtrophyPatterns.m
% Create a certain pattern:

clear all; clc;
addpath(genpath('/home/bkhanal/Documents/softwares/matlabTools/'));

% load b2d.mat;
% b3d = repmat(b2d, [1 1 10]);
% b3d = padarray(b3d,[0 0 1],0,'both');
% % imshow3D(b3d);
% save('b3d.mat','b3d');
% writemetaimagefile('bTest.mha',b3d,[1 1 1]);

load b3d.mat;
% writemetaimagefile('btestmha',b3d,[1 1 1]);
% imshow3D(b3d);

csfLeftVal = 50; csfCenterVal = 60; csfRightVal = 70; 

        bVal = 100;

csf = zeros(size(b3d));     b = csf;        a = csf;
csfLeft = csf;              csfCenter = csf; csfRight = csf; 

csfLeft(b3d==csfLeftVal) = 1;
csfCenter(b3d==csfCenterVal) = 1;
csfRight(b3d==csfRightVal) = 1;
csf(csfLeft | csfCenter | csfRight) = 1;

b(b3d == bVal) = 1;         %gray/white matter segmentation

% atrophy:
% Put atrophy everywhere in the brain, and compute it's sum.
a(b==1) = 0.2;
aInBrain = sum(a(:));

% CASE 1:
% Uniformly distribute it to the csf:
a1 = a;
a1(csf==1) = -aInBrain/sum(csf(:));
% Check that total atrophy should be zero:
display('a1 sum');
sum(a1(:))
writemetaimagefile('testAtrophy1.mha',a1,[1 1 1]);

% CASE 2:
% Zero Exapnsion on left:
a2 = a;
a2(csfCenter | csfRight) = -aInBrain/(sum(csfCenter(:)) + sum(csfRight(:)));
display('a2 sum:');
sum(a2(:))
writemetaimagefile('testAtrophy2.mha',a2,[1 1 1]);

% CASE 3:
% Zero Expansion on right:
a3 = a;
a3(csfLeft | csfCenter) = -aInBrain/(sum(csfLeft(:)) + sum(csfCenter(:)));
display('a3 sum:');
sum(a3(:))
writemetaimagefile('testAtrophy3.mha',a3,[1 1 1]);

% CASE 4:
% Zero Expansion on both sides:
a4 = a;
a4(csfCenter==1) = -aInBrain/sum(csfCenter(:));
display('a4 sum:');
sum(a4(:))
writemetaimagefile('testAtrophy4.mha',a4,[1 1 1]);
