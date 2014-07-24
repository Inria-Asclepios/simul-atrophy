% Create a test problem of size nXn.
% Will have a two isolated CSF region in the center, options for CSF to be
% connected or not connected.

function [seg_mask div_slice] = createSyntheticProblem2D(xnum,ynum,CONNECTED_CSF, maskType)
% create the Brain and Csf mask:
% Numbers of nodes
% xnum    =   45;             % Horizontal
% ynum    =   45;             % Vertical

% mask centers
m1CX = xnum/3;
m1CY = ynum/2;
m2CX = 2*m1CX;
m2CY = m1CY;

% mask radius
mR = 4;

% true values => Parenchyma, initialize with all true.
seg_mask = ~(false(xnum,ynum));

% Create CSF regions:
seg_mask(m1CY-mR:m1CY+mR, m1CX-mR:m1CX+mR) = false;
seg_mask(m1CY-mR:m1CY+mR, m2CX-mR:m2CX+mR) = false;

% Call a function:
% center_mask = [round(xnum/2) round(ynum/2)];
% brain_mask = circleMask(center_mask,xnum/4,xnum-1,ynum-1);
% brain_mask = ellipseMask(center_mask,7,4,xnum-1,ynum-1);

% divergence at the center of the cells
a = zeros(ynum,xnum);
a(int16(ynum/2)-1:int16(ynum/2)+1, int16(xnum/2)-1:int16(xnum/2)+1) = 0.2;
% a(int16(ynum/3)+2:int16(ynum/3)+5, int16(xnum/2)-1:int16(xnum/2)+1) = 0.2;
% a(2*int16(ynum/3):2*int16(ynum/3)+3, int16(xnum/2)-1:int16(xnum/2)+1) = 0.2;

div_slice = a;

% In the periphery add zeros of width two.
seg_mask = padarray(seg_mask,[5 5]);
div_slice = padarray(div_slice,[5 5]);
end
