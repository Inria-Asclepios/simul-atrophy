function [a_x a_y] = gradientAtEdges2D(a)

% Compute gradient which resides on the edges
a_x = a(:,2:end) - a(:,1:end-1);
a_x = padarray(a_x,[0 1], 'replicate', 'both');
a_y = a(2:end,:) - a(1:end-1,:);
a_y = padarray(a_y,[1 0], 'replicate', 'both');
end