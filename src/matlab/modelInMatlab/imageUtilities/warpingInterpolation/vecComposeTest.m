% vecComposeTest.ma = ones(5,5)*0.7;


clear all;
% ux = ones(5,5)*0.7;
% ux = padarray(ux,[2 2],0,'both');
% uy = ux;
% u(:,:,1)= ux;
% u(:,:,2) = uy;
% v = vecCompose(u,u);
% vv = expRecursive(u,3);


im = imread('cameraman.tif');
[m n] = size(im);

% Translation field
ux = ones(m-2,n-2)*30;
% ux = ones(m,n)*30;
ux = padarray(ux,[2 2],0,'both');
uy = ux;


% % Create rotation displacement field:
% theta_deg = 90;
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


u(:,:,1)= ux;
u(:,:,2) = uy;

%% 
uu = u;
im1 = im;
for indx = 1:1:10
    uu = expRecursive(indx*u,10);
    im1 = transformImage(im,uu);
    subplot(121), imshow(im,[]);
    subplot(122), imshow(im1,[]); title(['image ' num2str(indx)]);
    pause(0.1);
end
% im_out = transformImage(im,u);
% imshow(im), title('original image');
% figure, imshow(im_out, []), title('transformed image');

