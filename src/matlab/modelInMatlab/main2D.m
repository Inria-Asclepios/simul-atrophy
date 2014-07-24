%% Experiments
clear all; clc;
addpath(genpath('/home/bkhanal/works/AdLemModel/src/matlab'));
xnum = int32(45); ynum = int32(45); CONNECTED_CSF = false; mask_type = 'square';
[seg_mask div_slice] = createSyntheticProblem2D(xnum,ynum,CONNECTED_CSF,mask_type);

extra_CSF = 5; VENTRICLES_UNIFORM = true;
muBrain = 2500; muRatio = 10000;
lambdaBrain = 250; lambdaRatio = 10000;
[v p]=AdLemExperiments2D(seg_mask,div_slice,extra_CSF,muBrain,muRatio,lambdaBrain,lambdaRatio,VENTRICLES_UNIFORM);
% [v2 p2] = AdLemTaras2D(seg_mask,div_slice,true,muBrain,muRatio,lambdaBrain,lambdaRatio,true);

%% Load Image data:
clear all; clc;
addpath(genpath('/home/bkhanal/works/AdLemModel/src/matlab'));
dirname = '/home/bkhanal/works/AdLemModel/data/Template';

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
USE_REAL_DATA = 0;

%% Modify divergence data:
% 
% % % Diffusion
% % div_slice_b = diffuseAtrophy(div_slice,seg_mask,true);
% div_slice_csf = diffuseAtrophy(div_slice,~seg_mask,true);
% div_slice_m = div_slice;
% div_slice_m(~seg_mask) = 0;
% 

%% Get synthetic divmap:
% div_slice_m = div_slice;
div_slice_1 = getSynthAtrophy(div_slice,seg_mask,10);
subplot(121), imagesc(div_slice_1), axis image;
div_slice_2 = getSynthAtrophy(div_slice,seg_mask,400);
subplot(122), imagesc(div_slice_2), axis image;
%% Use the model to get velocity/displacement and pressure fields. 

muBrain = 2500;
lambdaBrain = 2000;

% muB/muc and lambdaB/lambdaC
muRatio = 1; lambdaRatio = 1;
% muRatio = 1000; lambdaRatio = 1000;
% muRatio = 2; lambdaRatio = 2;

% Augmented Lagrangian Method:
% [v2 p2] = AdLemAl2D(true,muBrain,muRatio,lambdaBrain,lambdaRatio,seg_mask,div_slice,USE_REAL_DATA);

% Coupled Method:
VENTRICLES_UNIFORM = true; extra_CSF = 2;
[v p]=AdLemExperiments2D(seg_mask,div_slice,extra_CSF,muBrain,muRatio,lambdaBrain,lambdaRatio,VENTRICLES_UNIFORM);
% [v2 p2] = AdLemTaras2D(seg_mask,div_slice_1,USE_REAL_DATA,muBrain,muRatio,lambdaBrain,lambdaRatio,true);
% [v3 p3] = AdLemTaras2D(seg_mask,div_slice_2,USE_REAL_DATA,muBrain,muRatio,lambdaBrain,lambdaRatio,true);

% [v2 p2] = AdLemTaras2D(seg_mask,div_slice_1,USE_REAL_DATA,muBrain,muRatio,lambdaBrain,lambdaRatio);
% [v3 p3] = AdLemTaras2D(seg_mask,div_slice_1,USE_REAL_DATA,muBrain,muRatio*1000000,lambdaBrain,lambdaRatio*1000000);


% imshow(seg_mask, []);

% saved_file = [dirname '/' cutting_plane num2str(slice_num) '_muR' num2str(muRatio) '.mat'];
% save(saved_file, 'img_slice', 'seg_mask','div_slice', 'v2', 'p2');
% staggeredAugLagr2D(true,muBrain,muRatio,lambdaBrain,lambdaRatio,seg_mask,div_slice,USE_REAL_DATA)
%% Load v2 and img_slice of certain muRatio:
clear v2; clear p;
muRatio = 1000; 
saved_file = [dirname '/' cutting_plane num2str(slice_num) '_muR' num2str(muRatio) '.mat'];
load(saved_file);

%% display the image with velocity superimposed:
figure,
% imshow(img_slice,[]);
vec_scale = 3000;
subplot(121), imshow(img_slice,[]);
hold on;
quiver(v2(:,:,1)*vec_scale,v2(:,:,2)*vec_scale,'Color','red','Autoscale','Off'); axis ij image; %use this instead of flipping!
title(['muBrain = ' num2str(muBrain) ',  muCSF = ' num2str(muBrain/muRatio)]);

vec_scale = 3000;
subplot(122), imshow(img_slice,[]);
hold on;
quiver(v3(:,:,1)*vec_scale,v3(:,:,2)*vec_scale,'Color','red','Autoscale','Off'); axis ij image; %use this instead of flipping!
title(['muBrain = ' num2str(muBrain) ',  muCSF = ' num2str(muBrain/(muRatio))]);
hold off;
figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',24,'fontWeight','bold')
figure,
subplot(121), imagesc(div_slice_1), axis image; colorbar;
subplot(122), imagesc(div_slice_2), axis image; colorbar;
figureHandle = gcf;
%# make all text in the figure to size 14 and bold
set(findall(figureHandle,'type','text'),'fontSize',24,'fontWeight','bold')

%% displaying 
img_slice1 = img_slice;
write_path = '/home/bkhanal/works/AdLemModel/results/fromMatlab/keBhayo';
% mov = avifile('tst_mov.avi');
figure,
sstp = 10;
num_it = 10;
for indx = sstp:sstp:num_it*sstp
%         v22 = expRecursive((1e3)*(indx-sstp)*v2,8);
        v22 = expRecursive((1e2)*(indx-sstp)*v2,8);
        v33 = expRecursive((1e2)*(indx-sstp)*v3,8);
    div_new2 = transformImage1(img_slice1,v22);
    imwrite(mat2gray(div_new2),[write_path 'Pointed_image' num2str(indx/sstp) '.png']);
%     imshow(div_new2,[]);
%     F = getframe(gca);
%     mov = addframe(mov,F);
    div_new3 = transformImage1(img_slice1,v33);
    imwrite(mat2gray(div_new3),[write_path 'Diffused_image' num2str(indx/sstp) '.png']);
%     subplot(121), imshow(img_slice,[]);
    subplot(121), imshow(div_new2,[]);
%     imshow(div_new3,[]);
    subplot(122), imshow(div_new3,[]), title(['image ' num2str(indx)]);
    pause(0.01);
%     pause(0.5);
%     img_slice1 = div_new;
end
% mov = close(mov);
%% Displaying four different time points:
figure, subplot(221), imshow(img_slice, []), title('baseline image');
img_slice1 = img_slice;
for indx_i = 1:3
    indx = 10*indx_i;
    v22 = expRecursive((1e3)*indx*v2,8);
    div_new = transformImage(img_slice1,v22);
    subplot(2,2,indx_i+1), imshow(img_slice1-div_new,[]), title(['image at ' num2str(indx) 't']);
end


%% Difference:
figure,
imshow(img_slice-div_new,[]), title('difference image');

%% video creation
mov = avifile('movie.avi');
count=0;
for i=1:10
    name1=strcat('video/image',num2str(i),'.png');
    a=imread(name1);
    while count<5
        count=count+1;
        gca = imshow(a);
        F=getframe(gca);
        mov=addframe(mov,F);
    end
    count=0;
end
close all
mov=close(mov);
