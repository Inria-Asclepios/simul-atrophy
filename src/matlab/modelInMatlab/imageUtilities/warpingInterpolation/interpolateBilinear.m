% Weighted bilinear interpolation:
% Input: 
% im_in --> input image
% row, col --> top-left discrete position: p
% dx, dy ---> distance of the point for which the interpolation is being
% done from p.
% --------------- p---------*------------------------------
% ----------------|   |     |------------------------------
% ----------------|   dy    |------------------------------
% ----------------|   |     |------------------------------
% ----------------|   x     |------------------------------
% ----------------|-dx-     |------------------------------
% ----------------*---------*------------------------------


function val_out = interpolateBilinear(im_in,row,col,dx,dy)

%outside four corners, so just copy corner value:
if (dx < 0 && dy < 0) 
    val_out = im_in(row,col);

% outside left or right vertical edges.
elseif (dx<0) 
    val_out = im_in(row,col)*(1-dy) + im_in(row+1,col)*dy;

% outside top or bottom horizontal edges
elseif (dy<0)
    val_out = im_in(row,col)*(1-dx) + im_in(row,col+1)*dx;

% Normal inner points
else
    val_out = im_in(row,col)*(1-dy)*(1-dx) + im_in(row+1,col)*dy*(1-dx)...
    + im_in(row,col+1)*(1-dy)*dx + im_in(row+1,col+1)*dy*dx;
end