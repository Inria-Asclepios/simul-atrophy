% Evolve segmented mask to have all the sources inside ventricles, i.e.
% effectively evolve the ventricles!


% Load Image data:
% clear all;
% addpath NIFTI_20110921/;
% dirname = 'RealData/Template';

% cutting_plane = 'coronal';
% slice_num = 116;
% [img_slice seg_mask div_slice] = getSegAndDivMask(dirname,cutting_plane,slice_num);

function div_slice_m = adaptDivMap(seg_mask, div_slice, num_it,dt)

div_slice = double(div_slice);
div_slice_mm = div_slice;
phi_o = ~seg_mask;

% Compute phi_o at grids from which gradients will be computed at center of
% the cells. Image has values at the center of the cells, so first need to
% compute the values at the edges:
disp_field = zeros(size(phi_o,1),size(phi_o,2),2);

figure;
for indx = dt:dt:num_it*dt+dt
%     Interpolate the values at staggered grid positions to comput gradients.
    phi_o_edge_ver = 0.5*(phi_o(:,1:end-1) + phi_o(:,2:end));
    phi_o_edge_hor = 0.5*(phi_o(1:end-1,:) + phi_o(2:end,:));
    % pad by replicating the border elements.
    phi_o_edge_hor = padarray(phi_o_edge_hor,[1 0],'replicate');
    phi_o_edge_ver = padarray(phi_o_edge_ver,[0 1],'replicate');

%     Compute gradient
    phi_x = phi_o_edge_ver(:,2:end) - phi_o_edge_ver(:,1:end-1);
    phi_y = phi_o_edge_hor(2:end,:) - phi_o_edge_hor(1:end-1,:);
    
    phi_grad_mag = sqrt((phi_x.^2) + (phi_y.^2));
        
    % Normalize it to obtain normal vector
    phi_nx = phi_x;
    phi_ny = phi_y;
    % Normalize all the non-zero vectors:
    nzero_indx = (phi_grad_mag ~= 0);
    phi_nx(nzero_indx) = phi_nx(nzero_indx)./(phi_grad_mag(nzero_indx));
    phi_ny(nzero_indx) = phi_ny(nzero_indx)./(phi_grad_mag(nzero_indx));
    
    phi_normal(:,:,1) = phi_nx;
    phi_normal(:,:,2) = phi_ny;

    %     weight them with div(v) values and time step:    
    phi_normal_w(:,:,1) = dt*(div_slice .* phi_normal(:,:,1));
    phi_normal_w(:,:,2) = dt*(div_slice .* phi_normal(:,:,2));
%     
    subplot(121), imagesc(phi_o), title('phi_o');
    hold on;
    quiver(phi_normal_w(:,:,1),phi_normal_w(:,:,2),'Color','yellow');
    axis ij image;
%     Evolve phi with weighted normal
    phi_t = transformImage1(phi_o,phi_normal_w);
%     Remask, expanding the CSF:
    phi_t(phi_t > 0) = 1;
    

    %     disp_field = vecCompose(invertDisplacementField(phi_normal_w),disp_field);
%     Compose the displacement field
    div_slice_mm = transformImageUinv(div_slice_mm,disp_field);
    disp_field = vecCompose(phi_normal_w,disp_field);
    
    subplot(122), imagesc(phi_t);
    hold on;
    quiver(phi_normal_w(:,:,1),phi_normal_w(:,:,2),'Color','yellow');
    axis ij image;
    title(['image iteration: ' num2str(indx/dt)]);

    pause(0.05);
    phi_o = phi_t;
end

div_slice_m = transformImageUinv(div_slice,disp_field);
clims = [min(min(div_slice(:)),min(div_slice_m(:))) max(max(div_slice(:)),max(div_slice_m(:)))];
figure,
subplot(121), imshow(div_slice,clims), title('original div map'), axis image;
subplot(122), imshow(div_slice_m,clims), title('div map modified'), axis image;

clims = [min(min(div_slice(:)),min(div_slice_mm(:))) max(max(div_slice(:)),max(div_slice_mm(:)))];
figure,
subplot(121), imshow(div_slice,clims), title('original div map'), axis image;
subplot(122), imshow(div_slice_mm,clims), title('div map modified m'), axis image;