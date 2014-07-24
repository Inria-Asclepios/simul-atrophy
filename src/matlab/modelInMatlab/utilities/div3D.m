% div3D

% For the end boundary of each dimension, zero is padded.
% A forward difference used.
% the last array element along each dimension of the result do not
% represent the true divergence unless the original input has zeros outside
% the border in reality too!!
% This is because, calculations at the last array elements along each
% dimensions are done assuming the input has zero everwhere outside the
% domain.

function div_a = div3D(a,h1,h2,h3)

[m n r d] = size(a);

% pad with zero after the last array element along each dimension
a_u = a(:,:,:,1);
a_u = padarray(a_u,[0 1 0],'post');

a_v = a(:,:,:,2);
a_v = padarray(a_v,[1 0 0],'post');

a_w = a(:,:,:,3);
a_w = padarray(a_w,[0 0 1],'post');

ux = a_u(:,2:n+1,:) - a_u(:,1:n,:);
vy = a_v(2:m+1,:,:) - a_v(1:m,:,:);
wz = a_w(:,:,2:r+1) - a_w(:,:,1:r);

div_a = (ux./h1) + (vy./h2) + (wz./h3);


end
