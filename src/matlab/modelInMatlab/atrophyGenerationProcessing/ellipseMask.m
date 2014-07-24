% Create a elliptical mask with logical true inside the ellipse and false
% outside the disk extending to a square of the size xnum X ynum.
% center represents center of the disk (x,y);

function mu = ellipseMask(center,radius1,radius2,xnum,ynum)

mu = logical(zeros(ynum,xnum));
for i = 1:ynum
    for j = 1:xnum
        if(((((i-center(2))^2)/(radius1^2) + ((j-center(1))^2)/(radius2^2)) - 1) <= 0)
            mu(i,j) = true;
        end
    end
end

end