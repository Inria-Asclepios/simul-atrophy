% gradient2D

% using forward difference
function [grad_a_x grad_a_y] = grad2D(a,h1,h2)

[m n] = size(a);

% dirichlet boundary condtion
a_ext = [a; 0.*a(m,:)];
a_ext = [a_ext 0.*a_ext(:,n)];

grad_a_x = a_ext(2:m+1,1:n)-a_ext(1:m,1:n); %IMP: vertical direction taken as x-axis!!!
grad_a_x = grad_a_x ./ h1;
%so that we can take 1st co-ordinate as x-axis when we write u(i,j)!!
grad_a_y = a_ext(1:m,2:n+1)-a_ext(1:m,1:n);
grad_a_y = grad_a_y ./ h2;

% neumann boundary condition
% a_ext = [a; a(m,:)];
% a_ext = [a_ext a_ext(:,n)];
% 
% grad_a_x = a_ext(2:m+1,1:n)-a_ext(1:m,1:n); %IMP: vertical direction taken as x-axis!!!
% grad_a_x = grad_a_x ./ h1;
% %so that we can take 1st co-ordinate as x-axis when we write u(i,j)!!
% grad_a_y = a_ext(1:m,2:n+1)-a_ext(1:m,1:n);
% grad_a_y = grad_a_y ./ h2;
end
