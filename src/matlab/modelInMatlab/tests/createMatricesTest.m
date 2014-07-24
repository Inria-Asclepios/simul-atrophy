% building different matrices that could be useful in solving partial
% differential equations.
% Tests!


m = 50;
n = 50;
h1 = 0.01;
h2 = 0.01;
N = 3*m*n;
% N = 2*m*n;
% f_tau = zeros(p,1);
f_tau = ones(N,1);

% Create the D2x and D2y matrices
l = [-2;1;zeros(n-2,1)];
g = [0;1;zeros(n-2,1)];
a = toeplitz(l,l);
b = toeplitz(-g,g);
D2x = 1/h1^2*kron(a,eye(m));
D2y = 1/h2^2*kron(eye(m),a);
Lp = D2x+D2y;


b = toeplitz(-g,g);
Dx = 1/h1*kron(eye(m),b);
Dy = 1/h2*kron(b,eye(m));

A = [Lp zeros(size(Lp)) Dx;
    zeros(size(Lp)) Lp Dy;
    Dx Dy zeros(size(Dx))];  %if -ve put then will be symmetric

% A = [Lp zeros(size(Lp));
%     zeros(size(Lp)) Lp];  %if -ve put then will be symmetric

% A_s = sparse(A);
x = A\f_tau;
% x = pcg(A_s,f_tau);
u = reshape(x(1:m*n),m,n);
v = reshape(x(m*n+1:2*m*n),m,n);
% p = reshape(x(2*m*n+1:3*m*n),m,n);
u = mat2gray(u);
imshow(u);
figure,
imagesc(A);
colormap(gray);

