% Compose two transormation by using displacement fields associated with
% the corresponding transformations.
% i.e. For T1oT2, inputs are v and w such that T1 = Id + v; T2 = Id + w.
% The output is u such that: T1oT2 = Id + u.
% 
% The deformation field corresponding to the composed transformation is:
% (VoW)(p) = W(p) + V(p + W(p)) for all points p in the image.
% 
% 
% V and W must be of the same size.
% V, W should have mXnX2 size with first and second matrices along 3rd
% dimension consisting of x-component and y-component of the velocity field
% respectively.
% 
%           o---->x-axis
%           |
%           |
%           v
%           y-axis
% 
% Bishesh Khanal,
% Asclepios, INRIA


function u = vecCompose(v,w)

vx = v(:,:,1); vy = v(:,:,2);
[rows cols] = size(vx);
wx = w(:,:,1); wy = w(:,:,2);

% Initialize the size of output
vwx = zeros(rows,cols);
vwy = vwx;

% Let p be a point at (i_indx,j_indx)
for i_indx = 1:rows
    for j_indx = 1:cols
%         if(i_indx == 5 && j_indx == 5)
%             test_debug = 1;
%         end
        %find p1 = p + w(p)
        %y-component is row direction
        p1_i = i_indx + wy(i_indx,j_indx); 

        if(p1_i > 1 && p1_i < rows)
            crow = floor(p1_i);
            dy = p1_i - crow;
        else
            dy = -1; %-ve an info to interpolation fn. 
            if(p1_i <= 1)
                crow = 1;
            else
                crow = rows;
            end
        end
        
        %x-component is col direction
        p1_j = j_indx + wx(i_indx,j_indx);
        if(p1_j > 1 && p1_j < cols)
            ccol = floor(p1_j);
            dx = p1_j - ccol;
        else
            dx = -1; %-ve an infor to interpolation fn.
            if(p1_j <= 1)
                ccol= 1;
            else
                ccol = cols;
            end
        end
        vwx(i_indx,j_indx) = wx(i_indx,j_indx) + interpolateBilinear(vx,crow,ccol,dx,dy);
        vwy(i_indx,j_indx) = wy(i_indx,j_indx) + interpolateBilinear(vy,crow,ccol,dx,dy);
        
    end
end

u(:,:,1) = vwx;
u(:,:,2) = vwy;
end
