% 
clear all;

n = 10;
img = zeros(n);
img(:,4) = 1;
% img(2,:) = 1;
uX = zeros(n);
uY = zeros(n);
% uy(4:end,1:n/2) = 2;
% uY(4,1:n/2) = 2;
% uX(1:n/2,4) = 0.086;
uX(1:n/2,4) = .7;
u(:,:,1) = uX;
u(:,:,2) = uY;

[img_n u_inv] = transformImage1(img,u);

img_n(img_n > 0 ) = 1;

figure,
subplot(121), imagesc(img);
hold on; quiver(uX,uY,'Color','yellow'), axis ij image;

ux = u_inv(:,:,1);
uy = u_inv(:,:,2);

subplot(122), imagesc(img_n);
hold on; quiver(ux,uy,'Color','yellow'), axis ij image;

% img = imread('cameraman.tif');
% [m n] = size(img);
% 
% % Translation field
% 
% % ux = ones(m,n)*30;
% % uy = ux;
% 
% 
% % % Create rotation displacement field:
% theta_deg = 30;
% rot_centre = [n/2 m/2];
% theta = (theta_deg/180)*pi;
% ux = zeros(m,n);
% uy = ux;
% for i_indx = 1:m
%     for j_indx = 1:n
%         ux(i_indx,j_indx) = (j_indx-rot_centre(1)) - ((j_indx-rot_centre(1))*cos(theta) - (i_indx-rot_centre(2))*sin(theta));
%         uy(i_indx,j_indx) = (i_indx-rot_centre(1)) - ((j_indx-rot_centre(1))*sin(theta) + (i_indx-rot_centre(2))*cos(theta));
%     end
% end
% 
% 
% u(:,:,1)= ux;
% u(:,:,2) = uy;
% 
% 
% % img_new = transformImage(img,u);
% img_new = transformImage1(img,u);
% 
% imshow(img,[]);
% figure, imshow(img_new,[]);
% 
% 
