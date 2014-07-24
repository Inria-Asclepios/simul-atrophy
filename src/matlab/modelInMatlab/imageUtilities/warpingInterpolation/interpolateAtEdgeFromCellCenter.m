% Interpolate values at the edge-centers from the values at cell centers.
% -----------------------------------------------------------------------
% -----------------------------------------------------------------------
% vertical downward => +ve y-axis, horizontal rightwards => +ve x-axis,
% horizontal towards the screen => +ve z-axis.
% -----------------------------------------------------------------------
% edge_axis = 1 or 2 or 3 corresponding to edge aligned to y, x or z-axis
% respectively.
% -----------------------------------------------------------------------
% If img is of size 4X4X4, for edge_axis = 2 (x-axis) we have the result
% at_edges of size 5X4X5 since the values for at_edges will be in stay in
% the yz-planes that cut the cube into half. Hence will have the same
% number of values as that if img along x-axis!!


function at_edges = interpolateAtEdgeFromCellCenter(img,edge_axis)
Y_AXIS = 1; X_AXIS = 2; Z_AXIS = 3;

% Pad the boundaries on all sides by copying the neighboring values.
img = cat(1,img(1,:,:),img,img(end,:,:));
img = cat(2,img(:,1,:),img,img(:,end,:));
img = cat(3,img(:,:,1),img,img(:,:,end));

% Interpolate
if(edge_axis == X_AXIS)
    at_edges = (img(:,:,1:end-1) + img(:,:,2:end))/2;
    at_edges = (at_edges(1:end-1,:,:) + at_edges(2:end,:,:))/2;
    %     Strip-off the paddings on x-axis:
    at_edges = at_edges(:,2:end-1,:);
elseif(edge_axis == Y_AXIS)
    at_edges = (img(:,:,1:end-1) + img(:,:,2:end))/2;
    at_edges = (at_edges(:,1:end-1,:) + at_edges(:,2:end,:))/2;
    %         Strip-off the paddings on y-axis:
    at_edges = at_edges(2:end-1,:,:);
elseif(edge_axis == Z_AXIS)
    at_edges = (img(1:end-1,:,:) + img(2:end,:,:))/2;
    at_edges = (at_edges(:,1:end-1,:) + at_edges(:,2:end,:))/2;
    %             Strip-off the paddings on z-axis:
    at_edges = at_edges(:,:,2:end-1);
else
    error('illegal argument for edge_axis');
end

end
