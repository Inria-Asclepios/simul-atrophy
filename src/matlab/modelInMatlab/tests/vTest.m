% vTest.m

clear vars;
% Generate a velocity field.


% grid size:
m = 20; n = 20;
c = [m/2-round(m/3); n/2-round(n/3)];   %center

% create the positions: horizontal - y-axis; vertical x-axis
[x, y] = meshgrid(0:n-1,0:m-1);


% create tralnslated positions
pos(:,:,1) = x-c(1);
pos(:,:,2) = y-c(2);  
% create distance function
% psi = 1./(1+sqrt(pos(:,:,1).^2 + pos(:,:,2).^2));
psi = 1;
% Now, compute velocity v:
sigma1 = -10;
sigma2 = -10;

M = [sigma1 0;
    0 sigma2];

u = psi .* (sigma1 .* x);
v = psi .* (sigma2 .* y);

u_t = psi .* (sigma1 .* pos(:,:,1));
v_t = psi .* (sigma2 .* pos(:,:,2));


% u_t = psi .* (sigma1 .* abs(pos(:,:,1)));
% v_t = psi .* (sigma2 .* abs(pos(:,:,2)));

% u_t = -psi*(sigma1 .* abs(pos(:,:,1)).*sign(pos(:,:,1)));
% v_t = -psi*(sigma2 .* abs(pos(:,:,2)).*sign(pos(:,:,2)));

% imagesc(psi);
% figure, quiver(u,v);


subplot(221), imagesc(pos(:,:,1)), title('x');
subplot(222), imagesc(pos(:,:,2)), title('y');
subplot(223), imagesc(u_t), title('u');
subplot(224), imagesc(v_t), title('v');

figure, quiver(u_t,v_t), axis ij image;



