% grad3D

% using forward difference (so as to make compatible with div3D for which
% backward difference was used, thus the laplacian is the centered one we
% get in taking div of grad!

% boundary condition we will take as dirichlet, i.e. zero outside the
% domain.

function [grad_x grad_y grad_z] = grad3D(a,h1,h2,h3)

[m n r d] = size(a);
% d not used but just to prevent r from having wrong value since a is 4D!
% as it is 3D image with 3 values at each voxels.


a_u = zeros(m+1,n,r);
a_u(1:end-1,:,:) = a;
grad_x = a_u(2:end,:,:) - a;
grad_x = grad_x ./ h1;
%IMP: vertical direction taken as x-axis!!!
%so that we can take 1st co-ordinate as x-axis when we write u(i,j,k)!!


a_v = zeros(m,n+1,r);
a_v(:,1:end-1,:) = a;
grad_y = a_v(:,2:end,:) - a;
grad_y = grad_y ./ h2;

a_w = zeros(m,n,r+1);
a_w(:,:,1:end-1) = a(:,:,:);
grad_z = a_w(:,:,2:end) - a;
grad_z = grad_z ./ h3;

end