% Find best transform:

function [img_new uxy] = transformImage1(img,u)
[rows cols] = size(img);
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

img_new = zeros(rows,cols);
        
for r_indx = 1:rows
    for c_indx = 1:cols
        %find X = x+ux
        %y-component is row direction
        X_i = r_indx + uy(r_indx,c_indx); 
        if(X_i > 1 && X_i < rows)
            cur_row = floor(X_i);
            dy = X_i - cur_row;
        else
            dy = -1; %-ve an info to interpolation fn. 
            if(X_i <=1)
                cur_row = 1;
            else
                cur_row = rows;
            end
        end
        
        %x-component is col direction
        X_j = c_indx + ux(r_indx,c_indx);
        if(X_j > 1 && X_j < cols)
            cur_col = floor(X_j);
            dx = X_j - cur_col;   

        else
            dx = -1; %-ve an info to interpolation fn.           
            if(X_j <=1)
                cur_col = 1;
            else
                cur_col = cols;
            end
        end
        img_new(r_indx,c_indx) = interpolateBilinear(img,cur_row,cur_col,dx,dy);
    end
end
        

end        