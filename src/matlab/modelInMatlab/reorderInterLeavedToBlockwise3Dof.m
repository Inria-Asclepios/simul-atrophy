% reorder the interleaved dof matrix to blockwise dof. i.e. if
% input is matrix acting on the vector [u1 v1 p1 u2 v2 p2 ... un vn pn] , 
% output will be a matrix acting on the vector [u1 u2 ... un v1 v2 ... vn
% p1 p2 ... pn].
% Currently works when each dof has equal number of entries, a square
% matrix input

function [resultA resultR] = reorderInterLeavedToBlockwise3Dof(A,R)
dof = 3;
u = A(1:dof:end,:);
v = A(2:dof:end,:);
p = A(3:dof:end,:);

uu = u(:,1:dof:end);
uv = u(:,2:dof:end);
up = u(:,3:dof:end);

vu = v(:,1:dof:end);
vv = v(:,2:dof:end);
vp = v(:,3:dof:end);

pu = p(:,1:dof:end);
pv = p(:,2:dof:end);
pp = p(:,3:dof:end);

resultA = [uu uv up;
    vu vv vp;
    pu pv pp];

resultR = [R(1:dof:end); R(2:dof:end); R(3:dof:end)];
end