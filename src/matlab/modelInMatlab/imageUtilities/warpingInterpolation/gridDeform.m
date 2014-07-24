% figure('Color', 'white')
% [X1_ndgrid,X2_ndgrid] = ndgrid(1:3,1:5);
% Z = zeros(3,5);
% mesh(X1_ndgrid,X2_ndgrid,Z)
% axis equal;
% view([0 0 1])

% clear all;

% ux= 0.5*ones(20);
% uy = 0.5*ones(20);
% ux(1:10,1:10) = 0;

% ux = v2(:,:,1);
% uy = v2(:,:,2);

% ux = upsample2D(ux,4,4);
% uy = upsample2D(uy,4,4);

% ux = ux*10;
% uy = uy*10;
% ux(1:30,:) = 0;
% uy(1:30,:) = 0;
% clear all;
% v = getRotField(4,4);

ux2 = v2(:,:,1)*2000;
uy2 = v2(:,:,2)*2000;

ux3 = v3(:,:,1)*2000;
uy3 = v3(:,:,2)*2000;


% ux_s = upsample2D(ux,3,3);
% uy_s = upsample2D(uy,3,3);

[ms ns] = size(ux3);
[x y] = meshgrid(1:ms,1:ns);
z = zeros(ms,ns);
% mesh(x,y,z), axis ij image;
% view([0 0 1]);

xx2 = x + ux2;
yy2 = y + uy2;

xx3 = x + ux3;
yy3 = y + uy3;

figure,
% subplot(121), quiver(ux,uy), axis ij image;
% subplot(122), mesh(xx,yy,z), axis ij image;

subplot(121), mesh(xx2,yy2,z), axis ij image; view([0 0 1]);
subplot(122), mesh(xx3,yy3,z), axis ij image; view([0 0 1]);



