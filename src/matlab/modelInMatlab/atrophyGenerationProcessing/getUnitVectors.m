% getUnitVectors of 3D vector field where vectors are arranged with thier
% components put in the first dimension. 2nd, 3rd and 4th dimension
% correspond to x, y and z axes correspondingly.

function out_field = getUnitVectors(in_field)

% Compute norm of each vectors.
rowNorm = sqrt(sum(in_field.^2,1));
% Divide every element of each vector with corresponding norms.
out_field = in_field ./ repmat(rowNorm,[3 1 1 1]);  
out_field(out_field==Inf) = 0;  
end