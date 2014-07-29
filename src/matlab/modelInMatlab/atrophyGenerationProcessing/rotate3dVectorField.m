% Rotate by given angles about the co-ordinate axes in the given order.
% Inputs:
% inp_vec_field : input vector field to be rotated.
% mask: TODO:(to be added), mask where the inp_vec_field is to be rotated.
% xAngle, yAngle, zAngle: rotation angles in radian about corresponding
% axes.
% order: order of the axes [x y z] about which the rotation will be
% applied. E.g. [2 1 3] will first rotate about y-axis followed by x-axis
% and finally z-axis.

function out_vec_field = rotate3dVectorField(inp_vec_field, xAngle, ...
    yAngle, zAngle, order)

inp_size = size(inp_vec_field);

Rx = [1 0 0
    0 cos(xAngle) -sin(xAngle)
    0 sin(xAngle) cos(xAngle)];

Ry = [cos(yAngle) 0 sin(yAngle)
    0 1 0
    -sin(yAngle) 0 cos(yAngle)];

Rz = [cos(zAngle) -sin(zAngle) 0
    sin(zAngle) cos(zAngle) 0
    0 0 1];

Rxyz(:,:,1) = Rx;
Rxyz(:,:,2) = Ry;
Rxyz(:,:,3) = Rz;

% rotation matrix:
R = eye(3);
for index = 1:3
    R = Rxyz(:,:,order==index)*R;
end

out_vec_field = zeros(size(inp_vec_field));

for i=1:inp_size(2)
    for j = 1:inp_size(3)
        for k = 1:inp_size(4)
            out_vec_field(:,i,j,k) = R * inp_vec_field(:,i,j,k);
        end
    end
end

end