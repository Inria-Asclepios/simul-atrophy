%% Create results: At ../../results/
clear all;clc;
addpath(genpath('/home/bkhanal/works/AdLemModel/src/matlab'));
res_path = '/user/bkhanal/home/works/AdLemModel/results/';

% fname = 'sol_cluster_big';
% sfName = 'size_cluster_big';

fname = 'sol';
sfName = 'size_lin_sys';

size_file = fopen([res_path sfName]); 
s = textscan(size_file,'%d',3);
s = s{1};
xm = s(1); ym = s(2); zm = s(3);
[ax ay az] = meshgrid(1:xm,1:ym,1:zm);

petscObj = PetscReadBinaryMatlab([res_path fname]);
vx = petscObj.x.vx;
vy = petscObj.x.vy;
vz = petscObj.x.vz;
p = petscObj.x.p;

vec3DToVtk(vx,vy,vz,ax,ay,az,[res_path 'velocity.vtk']);
savevtk(p,[res_path 'pressure.vtk']);

%% System matrix
clear all;clc;
addpath(genpath('/home/bkhanal/works/AdLemModel/src/matlab'));
res_path = '/user/bkhanal/home/works/AdLemModel/results/';

 sysFname = 'sys';
A = PetscBinaryRead([res_path sysFname]);

%% 
L = full(L);
A = full(A);
[rows cols] = size(A);
n = 4;
indxP = 1:(rows/n);
A_indx = zeros(1,rows);
A_indx(2:end) = 1:rows-1;
A_indx(:,(indxP-1)*4+1) = 4*indxP;
AA = A(A_indx,A_indx);

% length = m/4;
% Avx = zeros(length); Avy = Avx; Avz = Avy; Ap = Avz;
% indx = [0:length-1]';
% Avx = A(indx*4 + 1);