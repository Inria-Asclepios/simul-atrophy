% get two orthogonal vector fields for a given input vector field at each
% voxel.

function [v w] = getTwoOrthogonalVectorFields(u)
v = zeros(size(u));
% If performance is important check if mat2cell and arrayfun would be
% useful.
for i = 1:size(u,2)
    for j = 1:size(u,3)
        for k = 1:size(u,4)
            v(:,i,j,k) = getOrthogonalVector(u(:,i,j,k));
        end
    end
end
w = cross(u,v);
w = getUnitVectors(w);
end