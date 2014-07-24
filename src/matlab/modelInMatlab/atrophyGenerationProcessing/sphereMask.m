% Create a spherical mask with logical true on the sphere and false
% outside the sphere extending to a cube of the size xnum X ynum X znum.
% center represents center of the sphere: (x,y,z);

function mu = sphereMask(center,radius,xnum,ynum,znum)

mu = logical(zeros(ynum,xnum,znum));
for i = 1:ynum
    for j = 1:xnum
        for k = 1:znum
            if((((i-center(2))^2 + (j-center(1))^2 + (k-center(3))^2) - radius^2) <= 0)
                mu(i,j,k) = true;
            end
        end
    end
end

end