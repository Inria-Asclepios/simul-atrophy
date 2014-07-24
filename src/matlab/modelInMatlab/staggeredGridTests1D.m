% Staggered gridding nullspaces test 1D:
% Dirichlet Conditions:
clear all;clc;
vn = 3;
pn = 3;
A = toeplitz([1 zeros(1,vn)],[1 -2 1 zeros(1,vn-2)]);
A = A(1:vn,2:vn+1);

B = toeplitz([1 zeros(1,vn)],[1 -1 zeros(1,pn-1)]);
% B = B(1:vn,1:pn);
B = B(1:vn,vn-pn+1:end-(pn-vn+1));

S = [A B;
     B' zeros(pn)];
 rank(S)