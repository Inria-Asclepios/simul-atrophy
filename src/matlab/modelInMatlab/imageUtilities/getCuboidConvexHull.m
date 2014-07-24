% This function extracts a rectangular cuboid convex hull of logical true values
% in the given 3D images.
% Useful for example to extract rect. cuboid convex hull of brain parenchyma
% from a binary segmented images of brain parenchyma and other matters.
% ----------------------------------
% Inputs:
% inp_image :       an image with logical elements that has values either 0 or 1.
% ----------------------------------
% Outputs:
% img:              cuboidal convex hull.
% min_max_row_cols: array: [min_row max_row min_col max_col]
% --------------------------
% Bishesh Khanal
% Asclepios, INRIA Sophia Antipolis
% --------------------------

function [img min_max_y_x_z] = getCuboidConvexHull(inp_image)
% get the linear indices of all the positions that have value 1.
true_indx = find(inp_image);
[y x z] = ind2sub(size(inp_image),true_indx);
min_max_y_x_z = [min(y) max(y) min(x) max(x) min(z) max(z)];
img = inp_image(min_max_y_x_z(1) : min_max_y_x_z(2),...
    min_max_y_x_z(3) : min_max_y_x_z(4),...
    min_max_y_x_z(5) : min_max_y_x_z(6));
end