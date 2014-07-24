function u_inv = invertDisplacementField(u)
[rows cols] = size(u(:,:,1));
[X Y] = meshgrid(1:cols,1:rows);

uX = u(:,:,1);
uY = u(:,:,2);

phi_X = X + uX;
phi_Y = Y + uY;

ux_o = TriScatteredInterp(phi_X(:), phi_Y(:), -uX(:));
uy_o = TriScatteredInterp(phi_X(:), phi_Y(:), -uY(:));

ux = ux_o(X,Y);
uy = uy_o(X,Y);
uxy(:,:,1) = ux;
uxy(:,:,2) = uy;
u_inv = uxy;

end