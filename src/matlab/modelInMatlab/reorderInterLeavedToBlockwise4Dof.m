% reorder the interleaved dof matrix to blockwise dof. i.e. if
% input is matrix acting on the vector [u1 v1 p1 u2 v2 p2 ... un vn pn] , 
% output will be a matrix acting on the vector [u1 u2 ... un v1 v2 ... vn
% p1 p2 ... pn].
% Currently works when each dof has equal number of entries, a square
% matrix input

function [resultA resultR] = reorderInterLeavedToBlockwise4Dof(A,R)
dof = 4;
u = A(1:dof:end,:);
v = A(2:dof:end,:);
w = A(3:dof:end,:);
p = A(4:dof:end,:);

uu = u(:,1:dof:end);
uv = u(:,2:dof:end);
uw = u(:,3:dof:end);
up = u(:,4:dof:end);

vu = v(:,1:dof:end);
vv = v(:,2:dof:end);
vw = v(:,3:dof:end);
vp = v(:,4:dof:end);

wu = w(:,1:dof:end);
wv = w(:,2:dof:end);
ww = w(:,3:dof:end);
wp = w(:,4:dof:end);

pu = p(:,1:dof:end);
pv = p(:,2:dof:end);
pw = p(:,3:dof:end);
pp = p(:,4:dof:end);

resultA = [uu uv uw up;
    vu vv vw vp;
    wu wv ww wp;
    pu pv pw pp];

resultR = [R(1:dof:end); R(2:dof:end); R(3:dof:end); R(4:dof:end)];
end