% div2D
% axis:
% o-------> x (j)
% |
% |
% |
% y (i)

% using: (ax(i,j+1) - ax(i,j))/h1 + (ay(i+1,j) - ay(i,j))/h2 

% Note last row and col of the result contains values that may not be true
% since the boundaries of the input at the last row and col are extended by
% zeros before computing the divergence!


function div_a = div2D_new(a,h1,h2)

[m n d] = size(a);

% Assuming that at each position of 2D, the vectors are along the 3rd
% dimension, i.e. a(:,:,1) would represent 1st component of the vector

% add a row and col at the end with zeros
a_u = [a(:,:,1); 0.*(a(1,:,1))];
a_u = [a_u 0.*(a_u(:,1))];


%  add a row and col at the end with zeros.
a_v = [a(:,:,2); 0.*(a(1,:,2))];
a_v= [a_v 0.*(a_v(:,1))];


ux = a_u(1:m,2:n+1)-a_u(1:m,1:n); 
vy = a_v(2:m+1,1:n)-a_v(1:m,1:n);

div_a = (ux./h1) + (vy./h2);

end
