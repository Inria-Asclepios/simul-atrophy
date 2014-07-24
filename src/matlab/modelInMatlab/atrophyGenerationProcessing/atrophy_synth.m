% Synthetic atrophy:

%% Load Image data:
clear all;
addpath NIFTI_20110921/;
dirname = 'RealData/Template';

cutting_plane = 'coronal';
% slice_num = 110;
slice_num = 116;
% slice_num = 68;

% cutting_plane = 'axial';
% slice_num = 97;

% cutting_plane = 'sagittal';
% slice_num = ;


[img_slice seg_mask div_slice] = getSegAndDivMask(dirname,cutting_plane,slice_num);
% saved_file = [dirname '/' 'segAndDiv_' cutting_plane num2str(slice_num) '.mat'];
% save(saved_file, 'img_slice', 'seg_mask','div_slice');

% load(savedFile);
USE_REAL_DATA = 1;

%% Create synthetic atrophy
a = zeros(size(div_slice));


%Ventricles
a(86,78) = -10.6;
a(88,100) = -10.6;

%Hippocampus:
a(124,66) = 20.8;
a(124,67) = 20.8;
a(124,117) = 20.8;
a(124,118) = 20.8;

% Some other parts of brain:
a(61,69) = 30.8;
a(61,70) = 30.8;
a(60,117) = 30.8;
a(60,118) = 30.8;
figure,imagesc(a), axis ij image;

%% Blur the diracs:
% h = fspecial('gaussian',[20 20],5.8);
% a_b = imfilter(a,h);

% figure, imagesc(a_b), axis image;

%% Diffusion of a with diffusion coefficients Dc, Db (for CSF and Brain)
Dc = 1*(ones(size(seg_mask)).*(~seg_mask));
Db = 1*(ones(size(seg_mask)).*(seg_mask));

% Compute harmonic mean on edges:
Dcx = 2./((1./Dc(:,1:end-1)) + (1./Dc(:,2:end))); 
Dcy = 2./((1./Dc(1:end-1,:)) + (1./Dc(2:end,:))); 
% Extend them to edges by copying:
Dcx = padarray(Dcx,[0 1],'replicate','both');
Dcy = padarray(Dcy,[1 0],'replicate','both');

% Compute harmonic mean on edges:
Dbx = 2./((1./Db(:,1:end-1)) + (1./Db(:,2:end))); 
Dby = 2./((1./Db(1:end-1,:)) + (1./Db(2:end,:))); 
% Extend them to edges by copying:
Dbx = padarray(Dbx,[0 1],'replicate','both');
Dby = padarray(Dby,[1 0],'replicate','both');

%% Compute diffusion term:
a_diff = zeros(size(a));
dt = 0.2;


%% Brain region 
max_it = 500;
% subplot(121), imagesc(a), title('original image');
a_n = a;

for count = 1:max_it
    [a_x a_y] = gradientAtEdges2D(a_n);
    for ii = 1:size(a,1)
        for jj = 1:size(a,2)
            a_diff(ii,jj) = Dby(ii+1,jj)*a_y(ii+1,jj) - Dby(ii,jj)*a_y(ii,jj) ...
                + Dbx(ii,jj+1)*a_x(ii,jj+1) - Dbx(ii,jj)*a_x(ii,jj);
        end
    end
    a_n = a_n + dt*a_diff;
%     subplot(122), imagesc(a_n), title(['diffused image: ' num2str(count)]);
%     pause(2);
end

% close all;
% a_n(seg_mask) = a(seg_mask);
% figure, imagesc(a_n), axis ij image;


%% CSF region
max_it = 30;
for count = 1:max_it
    [a_x a_y] = gradientAtEdges2D(a_n);
    for ii = 1:size(a,1)
        for jj = 1:size(a,2)
            a_diff(ii,jj) = Dcy(ii+1,jj)*a_y(ii+1,jj) - Dcy(ii,jj)*a_y(ii,jj) ...
                + Dcx(ii,jj+1)*a_x(ii,jj+1) - Dcx(ii,jj)*a_x(ii,jj);
        end
    end
    a_n = a_n + dt*a_diff;
%     subplot(122), imagesc(a_n), title(['diffused image: ' num2str(count)]);
%     pause(2);
end

figure, imagesc(a_n), axis ij image;



