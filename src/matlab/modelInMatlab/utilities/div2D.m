% div2D
% axis:


% using backward difference (so as to make compatible with grad for which
% forward difference was used, thus the laplacian is the centered one we
% get in taking div of grad!

% boundary condition we will take as dirichlet, i.e. zero outside the
% domain.

function div_a = div2D(a,h1,h2)

[m n d] = size(a);

% Assuming that at each position of 2D, the vectors are along the 3rd
% dimension, i.e. a(:,:,1) would represent 1st component of the vector
a_u = [0.*(a(1,:,1)); a(:,:,1)];
a_u = [0.*(a_u(:,1)) a_u];


a_v = [0.*(a(1,:,2)); a(:,:,2)];
a_v = [0.*(a_v(:,1)) a_v];

ux = a_u(2:m+1,2:n+1)-a_u(1:m,2:n+1); %IMP: vertical direction taken as x-axis!!!
%so that we can take 1st co-ordinate as x-axis when we write u(i,j)!!

vy = a_v(2:m+1,2:n+1)-a_v(2:m+1,1:n);
div_a = (ux./h1) + (vy./h2);


% % if using neumann condition, i.e. derivative along the normal to the
% % boundary taken to be zero.
% a_u = [a(1,:,1); a(:,:,1)];
% a_u = [a_u(:,1) a_u];
% 
% 
% a_v = [a(1,:,2); a(:,:,2)];
% a_v = [a_v(:,1) a_v];
% 
% ux = a_u(2:m+1,2:n+1)-a_u(1:m,2:n+1); %IMP: vertical direction taken as x-axis!!!
% %so that we can take 1st co-ordinate as x-axis when we write u(i,j)!!
% 
% vy = a_v(2:m+1,2:n+1)-a_v(2:m+1,1:n);
% div_a = (ux./h1) + (vy./h2);
end
