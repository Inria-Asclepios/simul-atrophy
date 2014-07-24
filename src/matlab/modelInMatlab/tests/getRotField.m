% % Create rotation displacement field:

function u = getRotField(m,n)
% m = 20; n = 20;

theta_deg = 30;
rot_centre = [n/2 m/2];
theta = (theta_deg/180)*pi;
ux = zeros(m,n);
uy = ux;
for i_indx = 1:m
    for j_indx = 1:n
        ux(i_indx,j_indx) = (j_indx-rot_centre(1)) - ((j_indx-rot_centre(1))*cos(theta) - (i_indx-rot_centre(2))*sin(theta));
        uy(i_indx,j_indx) = (i_indx-rot_centre(1)) - ((j_indx-rot_centre(1))*sin(theta) + (i_indx-rot_centre(2))*cos(theta));
    end
end


u(:,:,1)= ux;
u(:,:,2) = uy;

% quiver(ux,uy), axis ij image;

end