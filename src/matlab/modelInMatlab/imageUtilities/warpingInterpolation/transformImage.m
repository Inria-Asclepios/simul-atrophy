
% my approach
function img_new = transformImage(img,u)
uX = u(:,:,1);
uY = u(:,:,2);

[Ynum Xnum] = size(img);
img_new = zeros(size(img));
wt = zeros(size(img));
for Y = 1:Ynum
    for X = 1:Xnum
        x = X + uX(Y,X);
        y = Y + uY(Y,X);
        x_int = floor(x);
        y_int = floor(y);

%         if inside the image space:
        if(((x_int > -1) && (x_int < (Xnum+1))) && ...
                ((y_int > -1) && (y_int < (Ynum+1))))
            dx = x - x_int;
            dy = y - y_int;
            pix_val = img(Y,X);
            if((x_int > 0 && y_int > 0) )%&& (dx < 0.5 && dy < 0.5))
                img_new(y_int,x_int) = img_new(y_int,x_int) + pix_val*(1-dx)*(1-dy);
                wt(y_int,x_int) = wt(y_int,x_int) + (1-dx)*(1-dy);                
            end
            if((x_int < Xnum && y_int > 0) )%&& (dx > 0.5 && dy < 0.5))
                img_new(y_int,x_int+1) = img_new(y_int,x_int+1) + pix_val*dx*(1-dy);
                wt(y_int,x_int+1) = wt(y_int,x_int+1) + dx*(1-dy);               
            end
            if((x_int > 0 && y_int < Ynum) )%&& (dx < 0.5 && dy > 0.5))
                img_new(y_int+1,x_int) = img_new(y_int+1,x_int) + pix_val*(1-dx)*dy;
                wt(y_int+1,x_int) = wt(y_int+1,x_int) + (1-dx)*dy;
            end
            if((x_int < Xnum && Ynum < Ynum) )%&& (dx > 0.5 && dy > 0.5))
                img_new(y_int+1,x_int+1) = img_new(y_int+1,x_int+1) + pix_val*dx*dy;
                wt(y_int+1,x_int+1) = wt(y_int+1,x_int+1) + dx*dy;            
            end
        end
            
    end
end

% normalize those which has weights.
% For those without any pixel coming from original, could be a good idea to
% preserve what original images had in those pixels or better could be to
% interpolate from the nearest neighbors of the new image.
nz_indx = find(wt);
img_new(nz_indx) = img_new(nz_indx) ./ wt(nz_indx);
% z_indx = 
end




            

        