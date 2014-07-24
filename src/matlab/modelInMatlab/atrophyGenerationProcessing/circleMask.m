% Create a disk mask with logical true on the disk and false
% outside the disk extending to a square of the size xnum X ynum.
% center represents center of the disk (x,y);

function mu = circleMask(center,radius,xnum,ynum)

mu = logical(zeros(ynum,xnum));
for i = 1:ynum
    for j = 1:xnum
        if((((i-center(2))^2 + (j-center(1))^2) - radius^2) <= 0)
            mu(i,j) = true;
        end
    end
end

end