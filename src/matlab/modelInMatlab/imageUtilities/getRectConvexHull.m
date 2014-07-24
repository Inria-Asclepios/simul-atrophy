% This function extracts a rectangular convex hull of logical true values
% in the given images.
% Useful for example to extract rectangular convex hull of brain parenchyma
% from a binary segmented images of brain parenchyma and other matters.
% ----------------------------------
% Inputs:
% inp_image :       an image with logical elements that has values either 0 or 1.
% ----------------------------------
% Outputs:
% img:              rectangular convex hull.
% min_max_row_cols: array: [min_row max_row min_col max_col]
% --------------------------
% Bishesh Khanal
% Asclepios, INRIA Sophia Antipolis
% --------------------------

function [img min_max_row_cols] = getRectConvexHull(inp_image)
% get the linear indices of all the positions that have value 1.
true_indx = find(inp_image);
[rows cols] = ind2sub(size(inp_image),true_indx);
min_max_row_cols = [min(rows) max(rows) min(cols) max(cols)];
img = inp_image(min_max_row_cols(1) : min_max_row_cols(2), min_max_row_cols(3) : min_max_row_cols(4));
end