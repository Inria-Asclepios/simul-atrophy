% restoreFromData3D.m

% View and recompute from the stored matlab data the results of
% linearElasticCSF3D.m

function restoreFromData3D(filename)

load(filename)

N = m*n*r;

u = zeros(m,n,r,3);
u(:,:,:,1) = reshape(x(1:N),m,n,r);
u(:,:,:,2) = reshape(x(N+1:2*N),m,n,r);
u(:,:,:,3) = reshape(x(2*N+1:3*N),m,n,r);
u1 = u(:,:,:,1);
u2 = u(:,:,:,2);
u3 = u(:,:,:,3);

p = reshape(x(3*N+1:4*N),m,n,r);
toc
% 
constraint = div3D(u,h1,h2,h3) + a;
slice_c = 4;
% slice_s = 4;
% slice_t = 4;
% 
subplot(221), imagesc(u1(:,:,slice_c)), title('u');
subplot(222), imagesc(u2(:,:,slice_c)), title('v');
subplot(223), imagesc(u3(:,:,slice_c)), title('w');
subplot(224), imagesc(p(:,:,slice_c)), title('p');
slice_c = 8;
figure,
subplot(221), imagesc(u1(:,:,slice_c)), title('u');
subplot(222), imagesc(u2(:,:,slice_c)), title('v');
subplot(223), imagesc(u3(:,:,slice_c)), title('w');
subplot(224), imagesc(p(:,:,slice_c)), title('p');


% [cay cax] = meshgrid(1:n,1:m); %coronal view
% [say sax] = meshgrid(1:r,1:m); %sagittal view
% [tay tax] = meshgrid(1:n,1:r); %top view, invert z-axis!!

[ax_y ax_x ax_z] = meshgrid(1:n,1:m,1:r);
vec3DToVtk(u1,u2,u3,ax_x,ax_y,ax_z,'letsee.vtk');
[ax_yb ax_xb ax_zb] = meshgrid(1:n-2*csf_w,1:m-2*csf_w,1:r-2*csf_w);
u1_b = u1(csf_w+1:m-csf_w,csf_w+1:n-csf_w,csf_w+1:n-csf_w);
u2_b = u2(csf_w+1:m-csf_w,csf_w+1:n-csf_w,csf_w+1:n-csf_w);
u3_b = u3(csf_w+1:m-csf_w,csf_w+1:n-csf_w,csf_w+1:n-csf_w);
vec3DToVtk(u1_b,u2_b,u3_b,ax_xb,ax_yb,ax_zb,'letseeBrain.vtk');
save_to_base(1);





