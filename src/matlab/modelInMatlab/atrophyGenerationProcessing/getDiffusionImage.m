% create tensor images from three eigenvectors and three eigenvalues:
% Inputs:
% three eigenvectors and correspondign eigenvalues.
% Size elements will be stored along the first dimension of the resulting
% image where elements are stored in lower triangular format i.e.
% [a11 a21 a22 a31 a32 a33]

function tensorImage = getDiffusionImage(u,v,w,lambda1,lambda2,lambda3)
tensorImage = zeros(6,size(u,2),size(u,3),size(u,4));

uSq = u.*u; 
vSq = v.*v; 
wSq = w.*w;

u1u2 = u(1,:,:,:) .* u(2,:,:,:,:);
v1v2 = v(1,:,:,:) .* v(2,:,:,:,:);
w1w2 = w(1,:,:,:) .* w(2,:,:,:,:);

u1u3 = u(1,:,:,:) .* u(3,:,:,:,:);
v1v3 = v(1,:,:,:) .* v(3,:,:,:,:);
w1w3 = w(1,:,:,:) .* w(3,:,:,:,:);

u2u3 = u(2,:,:,:) .* u(3,:,:,:,:);
v2v3 = v(2,:,:,:) .* v(3,:,:,:,:);
w2w3 = w(2,:,:,:) .* w(3,:,:,:,:);

tensorImage(1,:,:,:) = lambda1*uSq(1,:,:,:) + lambda2*vSq(1,:,:,:) + ...
    lambda3*wSq(1,:,:,:);

tensorImage(2,:,:,:) = lambda1*u1u2 + lambda2*v1v2 + lambda3*w1w2;

tensorImage(3,:,:,:) = lambda1*uSq(2,:,:,:) + lambda2*vSq(2,:,:,:) + ...
    lambda3*wSq(2,:,:,:);

tensorImage(4,:,:,:) = lambda1*u1u3 + lambda2*v1v3 + lambda3*w1w3;

tensorImage(5,:,:,:) = lambda1*u2u3 + lambda2*v2v3 + lambda3*w2w3;

tensorImage(6,:,:,:) = lambda1*uSq(3,:,:,:) + lambda2*vSq(3,:,:,:) + ...
    lambda3*wSq(3,:,:,:);

end
