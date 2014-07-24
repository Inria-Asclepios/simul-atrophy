% transformImage(img,u)
% Transforms a image with a given deformation field. 
% 1. Uses -u as the
% deformation field corresponding to the inverse transformation map.
% 
% 2. Uses bilinear interpolation to interpolate the values in img

function img_out = transformImageWrong(img,u)
ux = u(:,:,1);
uy = u(:,:,2);
[rows cols] = size(img);
img_out = zeros(rows,cols);
for r_indx = 1:rows
    for c_indx = 1:cols
        %find X = x-u
        %y-component is row direction
        X_i = r_indx - uy(r_indx,c_indx); 
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
        X_j = c_indx - ux(r_indx,c_indx);
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
        img_out(r_indx,c_indx) = interpolateBilinear(img,cur_row,cur_col,dx,dy);
    end
end

end