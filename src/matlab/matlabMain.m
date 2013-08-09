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
 sysFname = 'lin_sys';
A = PetscBinaryRead([res_path sysFname]);