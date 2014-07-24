% Apply diffusion or in simple case a laplace filter to the given region of
% interest in the input image. The region outside of interest will be
% set to zero. The function does not preserve the value of the input image
% outside the region of interest.
% Currently at the boundaries of region of
% interests, 
% i) symmetric ghost values is considered.
%
%
%
% 6th Feb, 2013
% Bishesh Khanal
% Asclepios, INRIA
%
function img_diffused = diffuseAtrophy(img, roi, DISPLAY)

% By default, do not display
if(nargin < 3)
    DISPLAY = false;
end


% temporary image to work with, set zero values outside roi, and extend it
% by padding it symmetrically with zeros of width 1.
img_ext = img;
img_ext(~roi) = 0;

% pad the image all around symmetrically with zeros, so that we don't
% exceed the index.
img_ext  = padarray(img_ext,[1 1]);
% extended roi
roi_ext = padarray(roi,[1 1]);


% Get the indices of roi_ext.
[row col] = find(roi_ext);

% inverted mask of the region of interest
roi_inv_ext = ~roi_ext;

laplacian= zeros(size(img_ext));

num_iter = 2;
for iter_indx = 1:num_iter
    for indx = 1:1:size(row)
        ci = row(indx); cj = col(indx);
        % Apply the Laplacian filter (homogenous isotropic diffusion):
        cv = img_ext(ci,cj); %current value being assigned to out-of-region pixels.
        laplacian(ci,cj) = (img_ext(ci-1,cj) + cv*roi_inv_ext(ci-1,cj))...
            + (img_ext(ci+1,cj) + cv*roi_inv_ext(ci+1,cj))...
            + (img_ext(ci,cj-1) + cv*roi_inv_ext(ci,cj-1))...
            + (img_ext(ci,cj+1) + cv*roi_inv_ext(ci,cj+1))...
            - 4*cv;
    end
    % Probably there is a need of normalization here!!
    img_ext = img_ext + laplacian;
    if(DISPLAY == true)
        subplot(121), imagesc(1000*img);
        axis image;
        subplot(122), imagesc(1000*img_ext);
        axis image;
%         pause(2);
    end
end
img_diffused = img_ext(2:end-1,2:end-1);

end