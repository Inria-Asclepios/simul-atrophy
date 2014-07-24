% Adaptive Augmented Lagrangian method for high discontinuity in viscosity.
% Solution of the Brain deformation with irregular domain (using extended
% domain method i.e. brain-CSF considered as an interface and extending the
% boundary to the skull.
% Regular staggered grid using a pressure-velocity formulation
% for a medium with varying lame's parameters.
%
function [L p] = AdLemAl2D(BUILD_MATRIX,muBrain,muRatio,lambdaBrain,lambdaRatio,seg_mask,div_slice,USE_REAL_DATA)

if(nargin < 8)
    USE_REAL_DATA = false;
    
end


if (USE_REAL_DATA == true)
    % Since mu and lambda Brain are greater than CSF, the ratio muB/muC given
    % as input in the form muRatio should be greater or equal to 1.
    if(lambdaRatio < 1 || muRatio < 1)
        error('invalid muRatio or invalid lamdaRatio, they should be greater than one');
    end
    
    [brain_mask min_max_row_col] = getRectConvexHull(seg_mask);
    [yn xn] = size(brain_mask);
    
    a = div_slice(min_max_row_col(1) : min_max_row_col(2), min_max_row_col(3) : min_max_row_col(4));
    %pad with extra_CSF width of zeros all around.
    extra_CSF = 2;
    brain_mask = padarray(brain_mask,[extra_CSF extra_CSF]);
    
    % distributed source outside convex hull(padded areas), excluding 4
    % corners where CC is not enforced! Total pixels of padded areas of
    % width w is (xn+2w)(yn+2w)-xn.yn
    a = padarray(a,[extra_CSF extra_CSF],-sum(a(:))/(2*extra_CSF*(xn + yn + 2*extra_CSF)));
    
    %Grid nodes have one dimension greater than the cell centers!
    yn = yn+extra_CSF*2+1;
    xn = xn+extra_CSF*2+1;
    
    
    % add a dummy atrophy corresponding to dummy pressure values!
    a = [zeros(1,xn); zeros(yn-1,1) a];
    
    %     brain mask ma dummy add garera vx jastai banaune, ani pheri add garne
    %     crop gareka kura haru, tyaspachhi seg_mask ra yo identical hunchha ki
    %     hunna check garne!!!
    %     brain_mask_test = [zeros(1,xnum); zeros(ynum-1,1) brain_mask];
    %     [m n] = size(div_slice);
    %     vc = size(brain_mask_test,2);
    %     brain_mask_test = cat(1,zeros(min_max_row_col(1)-extra_CSF-1-1, vc),brain_mask_test,zeros(m-min_max_row_col(2)-extra_CSF-1+1,vc));
    %     brain_mask_test = cat(2,zeros(m, min_max_row_col(3)-extra_CSF-1-1),brain_mask_test,zeros(m, n-min_max_row_col(4)-extra_CSF-1+1));
    
    
    
    
    % Model size from voxel resolution : xres X yres X zres
    % For 1mm X 1mm X 1mm
    xres = 1e-3; yres = 1e-3;
    xsize   =   (xn-1)*xres;        % Horizontal
    ysize   =   (yn-1)*yres;        % Vertical
    
else
    % Model size, m
    xsize   =   1;        % Horizontal
    ysize   =   1;        % Vertical
    
    % Numbers of nodes
    xn    =   21;             % Horizontal
    yn    =   21;             % Vertical
    % Call a function:
    center_mask = [round(xn/2) round(yn/2)];
    % brain_mask = circleMask(center_mask,xnum/4,xnum-1,ynum-1);
    % brain_mask = ellipseMask(center_mask,xnum/3,xnum/4,xnum-1,ynum-1);
    % brain_mask = ellipseMask(center_mask,7,4,xnum-1,ynum-1);
    brain_mask = circleMask(center_mask,xn/3,xn-1,yn-1);
    % brain_mask = circleMask(center_mask,7,xnum-1,ynum-1);
    
    % atrophy at the center of the cells
    a = computeAtrophy(yn-1,xn-1);
%     a(brain_mask==false) = 0;
    % distributed source outside brain,
    a(brain_mask==false) = -sum(a(:))/((xn-1)*(yn-1) - sum(brain_mask(:)));
%     a(1,1) = -sum(a(:));
    % add a dummy atrophy corresponding to dummy pressure values!
    a = [zeros(1,xn); zeros(yn-1,1) a];
    
    
end

% Grid step
xstp    =   xsize/(xn-1); % Horizontal
ystp    =   ysize/(yn-1); % Vertical

% Create vectors for nodal points positions (basic nodes)
% x       =   0:xstp:xsize;   % Horizontal
% y       =   0:ystp:ysize;   % Vertical

% Create vectors for cell centers positions (staggered nodes)
% xc      =   xstp/2:xstp:xsize-xstp/2; % Horizontal
% yc      =   ystp/2:ystp:ysize-ystp/2; % Vertical


% Model viscosity
muCSF = muBrain/muRatio;
lambdaCSF = lambdaBrain/lambdaRatio;

% mu and lambda in center of the cells!
muc = zeros(yn-1,xn-1);
lambdac = zeros(yn-1,xn-1);

muc(brain_mask==true) = muBrain;
muc(brain_mask==false) = muCSF;
lambdac(brain_mask==true) = lambdaBrain;
lambdac(brain_mask==false) = lambdaCSF;

% Interpolate mu and lambda on basic nodes from values at centers of the cells

mug = (muc(:,1:xn-2) + muc(:,2:xn-1))/2; %horizontal interpolation
mug = (mug(1:yn-2,:) + mug(2:yn-1,:))/2; %vertical interpolation

% Pad the values on the boundary edges with that of adjacent ones.
% horizontal padding
mug = [mug(1,:); mug; mug(end,:)];
% vertical padding
mug = [mug(:,1) mug mug(:,end)];


% mu at vertical edges: muVer, and horizontal muHor (useful for creating
% RHS)
% muVer = (muc(:,1:xn-2) + muc(:,2:xn-1))/2; %horizontal interpolation
% Pad all around by copying the adjacent values:
% muVer = padarray(muVer,[1 1],'replicate');

% muHor = (muc(1:yn-2,:) + muc(2:yn-1,:))/2; %vertical interpolation
% muHor = padarray(muHor,[1 1],'replicate');


% Now Add dummy zeros on 1st row and 1st col of centered cells:
muc = [zeros(1,xn); zeros(yn-1,1) muc];



lambdag = (lambdac(:,1:xn-2) + lambdac(:,2:xn-1))/2; %horizontal interpolation
lambdag = (lambdag(1:yn-2,:) + lambdag(2:yn-1,:))/2; %vertical interpolation

% Pad the values on the boundary edges with that of adjacent ones.
% horizontal padding
lambdag = [lambdag(1,:); lambdag; lambdag(end,:)];
% vertical padding
lambdag = [lambdag(:,1) lambdag lambdag(:,end)];

% Now Add dummy zeros on 1st row and 1st col of centered cells:
lambdac = [zeros(1,xn); zeros(yn-1,1) lambdac];

% Adaptive Lagrange parameter:
rc = ones(yn-1,xn-1);
% If not Adaptive: rc = 180*rc;
rc(brain_mask==true) = muBrain/10;
rc(brain_mask==false) = muCSF*muRatio/10;

% rc(brain_mask==true) = 1/muBrain;
% rc(brain_mask==false) = 1/muCSF;

% Compute at the edges:
ry = (rc(1:end-1,:) + rc(2:end,:))/2;
ry = padarray(ry,[1 0],'replicate','symmetric');
ry = padarray(ry,[0 1],'replicate','post');

rx = (rc(:,1:end-1) + rc(:,2:end))/2;
rx = padarray(rx,[0 1],'replicate','symmetric');
rx = padarray(rx,[1 0],'replicate','post');

% Now extend a dummy on the 1st row and 1st col of centered cells:
rc = padarray(rc,[1 1], 'replicate','pre');

% Computing Kcont and Kbond coefficients
kcont   =   2*muBrain/(xstp+ystp);
kbond   =   4*muBrain/(xstp+ystp)^2;
% kcont   =   2*muCSF/(xstp+ystp);
% kbond   =   4*muCSF/(xstp+ystp)^2;

DIRICHLET_BC = 0;
NEUMANN_BC = 1;
PSI_BC     = 2;

RELAX_CC = 1;
ENFORCE_CC = 0;

% Boundary Condition Options:
bc = DIRICHLET_BC; %dirichlet condition on displacement.
% bc = NEUMANN_BC; %neumann condition on displacement.
% bc = PSI_BC;   %psi function boundary condition on displacement.

% Pressure condition in one cell (i,j)
% p0cell  =   20;

% Options for enforcing pressure value at a point
% 1st option: Relax the CC at that point and fix the pressure value, this
% way can be done by simply replacing the row corresponding to CC
% cosnstraint by a row corresponding to setting specific value for
% pressure.
% 2nd option: Do not change the CC, instead add extra row specifying the
% pressure value at the corresponding point.

pChoice = RELAX_CC;  %By relaxing compressibility condition at p0cell, eqv.
% done in the code by replacing corresponding CC row by set pressure at corresponding point.

% pChoice = ENFORCE_CC;  %By enforcing compressibility condition at p0cell,
% means we need to add extra row to fix the pressure at that point.

% Distributed sources:
% pSourcesD = true;   %Create distributed sources by relaxing CC at certain pts.
% pSourcesD = false;   %No distributed sources.

% Number of Extra constraints required, 2 for translation, 1 for rotation
% if (bc == NEUMANN_BC)
if (bc == NEUMANN_BC)
    ecnstrnt = 3;
else
    ecnstrnt = 0;
end

if(pChoice == ENFORCE_CC)   %if CC enforced, extra constraints on pressure required to make full rank matrix.
    ecnstrnt = ecnstrnt + 1;
end

% if(pSourcesD == true )   %If distributed sources created by relaxing CC,
% %     pressure must be specified, which is of course done inside the loop by row replacing, so this addition of the constraint should not be needed I guess!
%     ecnstrnt = ecnstrnt + 4;
% end
% index keeping track of the ecnstrnt
% ecnstrnt_indx = 1;


% Top left:
% p0cellx = 3; %indx j
% p0celly = 2; %indx i

% Top right:
% p0cellx = xnum-1; %indx j
% p0celly = 2; %indx i

% Bottom left:
% p0cellx = 3; %indx j
% p0celly = ynum-1; %indx i

% Center
% p0cellx = round(xnum/2); %indx j
% p0celly = round(ynum/2); %indx i

% Bottom right
% p0cellx = xnum-1; %indx j
% p0celly = ynum-1; %indx i

% Arbitrary point:
% p0cellx = round(xnum/2)+4;
% p0celly = round(ynum/2)-4;


% Get the atrophy data at the cell-centers

if(bc == PSI_BC)
    % Solve for psi s.t. v = w - grad(psi), i.e. solve for:
    % psi_xx + psi_yy = a;
    psi = solvePsi(a,xstp,ystp);
    % strip off zeros, and make psi of the same size
    psi = psi(2:yn+1,2:xn+1);
end


% Matrix of coefficients initialization
% L       =   sparse(ecnstrnt+xnum*ynum*3,xnum*ynum*3);
% Tentative number of non-zeros in L matrix
% Needs to change later, slightly bigger right now:
tot_non_zeros = 3*(xn+yn) + 18*((yn-3)*(xn-2) + ...
    (xn-3)*(yn-2)) + 5 + 2*yn + 2*xn; %last ko 2*yn+2*xn ettikai extra
i_indx = zeros(tot_non_zeros,1);
j_indx = zeros(tot_non_zeros,1);
L_vals = zeros(tot_non_zeros,1);


% Vector of right part initialization
% R       =   zeros(ecnstrnt+ xn*yn*2,1);


% Solving x-Stokes, y-Stokes and continuity equations
% x-Stokes: mu(d2vx/dx2+d2vx/dy2)-dP/dx=(mu+lambda)da/dx
% y-Stokes: mu(d2vy/dx2+d2vy/dy2)-dP/dy=gy*RHO + (mu+lambda)da/dy
% continuity: dvx/dx+dvy/dy=-a
% Compose matrix of coefficients L()
% and vector (column) of right parts R()
% Boundary conditions: free slip
% Process all Grid points


vals_cntr = 0;
    function Lf(ii,jj,xx)
        vals_cntr = vals_cntr + 1;
        i_indx(vals_cntr) = ii;
        j_indx(vals_cntr) = jj;
        L_vals(vals_cntr) = xx;
    end

%
% initialize p as zero with dummy, but when computing p_n+1 from p_n, use
% without dummy.
boundary_mask = false(yn,xn);
boundary_val = zeros(yn,xn);

    function R = setRhs(p_in)
        R       =   zeros(ecnstrnt+ xn*yn*2,1);
        for ii = 1:1:yn
            for jj = 1:1:xn
                iinvx    =   ((jj-1)*yn+ii)*2-1; % invx
                iinvy    =   iinvx+1;
                % Ghost vx unknowns (i=ynum) and boundary nodes (i=1, i=ynum-1, j=1, j=xnum)
                if(ii==1 || ii==yn-1 || ii==yn || jj==1 || jj==xn)
                    if(ii==yn)
                        R(iinvx,1) = 0;
                    end
                    if((jj==1 ) && ii<yn)
                        if (bc == DIRICHLET_BC)
                            R(iinvx,1)   = 0;
                        elseif (bc == NEUMANN_BC)
                            % Neumann condition: -vx(i,j+1)+vx(i,j)=0
                            R(iinvx,1)           =   0;          % Right-hand-side part
                            %                     elseif (bc == PSI_BC)
                            %psi boudnary condition
                            %                         R(iinvx,1) = kbond*(psi(ii,jj)/xstp);
                        end
                    end
                    %             Right boundary
                    if((jj==xn) && ii<yn)
                        if(bc==DIRICHLET_BC)
                            %                     Dirichlet condition: vx(i,j) = 0
                            R(iinvx,1)    = 0;
                        elseif(bc==NEUMANN_BC)
                            % %                 Neumann condition: vx(i,j) - vx(i,j-1) = 0
                            R(iinvx,1) = 0;
                            %                 elseif(bc==PSI_BC)
                            %             psi boundary condition
                            %                     R(iinvx,1) = kbond*(psi(i,j)/xstp);
                        end
                    end
                    % Upper boundary, inner points (i=1, 1<j<xnum)
                    if(ii==1 && jj>1 && jj<xn)
                        if(bc==DIRICHLET_BC)
                            %                     Dirichlet condition: vx(i,j) = 0
                            % % equivalent to No slip vx=0: vx(i,j)-1/3*vx(i+1,j)=0.
                            R(iinvx,1)          =   0;          % Right-hand-side part
                        elseif(bc==NEUMANN_BC)% || bc == PSI_BC)
                            % Neumann condition dvx/dy=0: vx(i+1,j)-vx(i,j)=0
                            %                 Same for psi boundary condition!
                            R(iinvx,1)           =   0;          % Right-hand-side part
                        end
                    end
                    % Lower boundary, inner points (i=ynum-1, 1<j<xnum)
                    if(ii==yn-1 && jj>1 && jj<xn)
                        if(bc==DIRICHLET_BC)
                            %                     Dirichlet or No slip:
                            % % No slip vx=0: vx(i,j)-1/3*vx(i-1,j)=0
                            R(iinvx,1)         =   0; % Right part
                        elseif(bc==NEUMANN_BC)% || bc == PSI_BC)
                            % Neumann condition dvx/dy=0: vx(i,j)-vx(i-1,j)=0
                            R(iinvx,1)           =   0; % Right-hand-side part
                        end
                    end
                    %Internal nodes: mu(d2vx/dx2+d2vx/dy2) + r d/dx(dvx/dx +
                    %dvy/dy) = dP/dx + (r + mu + lambda)da/dx
                else
                    % Right part:dp_in/dx + (r+mu+lambda)da/dx =
                    % (p(i+1,j+1) - p(i+1,j))/dx + const.*(a(i+1,j+1)-a(i+1,j))/dx
                    R(iinvx,1) = (p_in(ii+1,jj+1) - p_in(ii+1,jj))/xstp + ...
                        (-rx(ii,jj) + (muc(ii+1,jj+1)+muc(ii+1,jj))/2 + (lambdac(ii+1,jj+1)+ ...
                        lambdac(ii+1,jj))/2) * (a(ii+1,jj+1)-a(ii+1,jj))/xstp;
                end
                
                % y-Stokes equation
                % Ghost vy unknowns (j=xnum) and boundary nodes (i=1, i=ynum, j=1, j=xnum-1)
                if(ii==1 || ii==yn || jj==1 || jj==xn-1 || jj==xn)
                    % Ghost vy unknowns (j=xnum: vy(i,j)=0
                    if(jj==xn)
                        R(iinvy,1)           =   0;
                    end
                    % Upper boundary, i=1
                    if((ii==1) && (jj<xn))
                        if(bc==DIRICHLET_BC)
                            %                     Dirichlet condition vy = 0
                            R(iinvy,1) = 0;
                        elseif(bc==NEUMANN_BC)
                            % %            Neumann codition: -vy(i+1,j) + vy(i,j) = 0;
                            R(iinvy,1)           =   0;
                            %                 elseif(bc==PSI_BC)
                            %             psi boundary condtion:
                            %                     Lf(invy,invy,1*kbond);
                            %                     R(invy,1) = kbond*(psi(i,j)/ystp);
                        end
                    end
                    %             Lower boundary, i=ynum
                    if(ii==yn && jj < xn)
                        if(bc==DIRICHLET_BC)
                            %                     Dirichlet condition vy = 0
                            R(iinvy,1) = 0;
                        elseif(bc==NEUMANN_BC)
                            % %           Neumann codition: vy(i,j) - vy(i-1,j) = 0;
                            R(iinvy,1) = 0;               %Right side
                            %                 elseif(bc==PSI_BC)
                            %                 psi boundary condtion:
                            %                     Lf(invy,invy,1*kbond);
                            %                     R(invy,1) = kbond*(psi(i,j)/ystp);
                        end
                    end
                    
                    % Left boundary, inner points (j=1, 1<i<ynum)
                    if(jj==1 && ii>1 && ii<yn)
                        if(bc==DIRICHLET_BC)
                            %                     Dirichlet or eqv. to no slip:
                            %             % No slip vy=0: vy(i,j)-1/3*vy(i,j+1)=0
                            R(iinvy,1)=0;
                        elseif(bc==NEUMANN_BC)% || bc == PSI_BC)
                            % Neumann dvy/dx=0: vy(i,j)-vy(i,j+1)=0
                            %                 eqv. for psi boundary condtion
                            R(iinvy,1)=0;
                        end
                    end
                    % Right boundary, inner points (j=xnum-1, 1<i<ynum)
                    if(jj==xn-1 && ii>1 && ii<yn)
                        if(bc==DIRICHLET_BC)
                            %                 Dirichlet or no slip: vy=0: vy(i,j)-1/3*vy(i,j-1)=0
                            R(iinvy,1)=0;
                        elseif(bc==NEUMANN_BC )%|| bc==PSI_BC)
                            % Free slip dvy/dx=0: vy(i,j)-vy(i,j-1)=0
                            %                 eqv. to psi boundary condition
                            R(iinvy,1)=0;
                        end
                    end
                    %Internal nodes
                else
                    % Right part: dp/dy + (r+mu+lambda)da/dy
                    R(iinvy,1) = (p_in(ii+1,jj+1) - p_in(ii,jj+1))/ystp + ...
                        (-ry(ii,jj) + (muc(ii+1,jj+1)+muc(ii,jj+1))/2 + (lambdac(ii+1,jj+1) + ...
                        lambdac(ii,jj+1))/2)*(a(ii+1,jj+1)-a(ii,jj+1))/ystp;
                end
            end
        end
    end


if (BUILD_MATRIX == true)
tic;
for i=1:1:yn
    for j=1:1:xn
        
        % Global index for vx, vy in the current node
        invx    =   ((j-1)*yn+i)*2-1; % invx
        invy    =   invx+1;
        
        % x-Stokes equation
        % Ghost vx unknowns (i=ynum) and boundary nodes (i=1, i=ynum-1, j=1, j=xnum)
        if(i==1 || i==yn-1 || i==yn || j==1 || j==xn)
            % Ghost vx unknowns (i=ynum: vx(i,j)=0
            if(i==yn)
                Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                %                 R(invx,1   )        =   0;          % Right-hand-side part
                boundary_mask(i,j) = true;
            end
            % Left and Right boundaries (j=1, j=xnum)
            %             Left boundary
            if((j==1 ) && i<yn)
                if (bc == DIRICHLET_BC)
                    %                 Dirichlet condition: vx(i,j) = 0
                    Lf(invx,invx,1*kbond);  %Coeff. for vx(i,j)
                    %                     R(invx,1)   = 0;
                    boundary_mask(i,j) = true;
                elseif (bc == NEUMANN_BC)
                    % Neumann condition: -vx(i,j+1)+vx(i,j)=0
                    Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                    Lf(invx,invx+yn*2,-1*kbond);    % coeff. for vx(i,j+1)
                    %                     R(invx,1)           =   0;          % Right-hand-side part
                    boundary_mask(i,j) = true;
                    %
                elseif (bc == PSI_BC)
                    %                 psi boudnary condition
                    Lf(invx,invx,1*kbond);
                    %                     R(invx,1) = kbond*(psi(i,j)/xstp);
                    boundary_mask(i,j) = true;
                    boundary_val(i,j) = kbond*(psi(i,j)/xstp);
                end
            end
            %             Right boundary
            if((j==xn) && i<yn)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition: vx(i,j) = 0
                    Lf(invx,invx,1*kbond); %Coeff. for vx(i,j)
                    %                     R(invx,1)    = 0;
                    boundary_mask(i,j) = true;
                elseif(bc==NEUMANN_BC)
                    % %                 Neumann condition: vx(i,j) - vx(i,j-1) = 0
                    Lf(invx,invx,1*kbond);  %Coeff. for vx(i,j)
                    Lf(invx,invx-yn*2,-1*kbond); %Coeff. for vx(i,j-1)
                    %                     R(invx,1) = 0;
                    boundary_mask(i,j) = true;
                elseif(bc==PSI_BC)
                    %             psi boundary condition
                    Lf(invx,invx,1*kbond);
                    %                     R(invx,1) = kbond*(psi(i,j)/xstp);
                    boundary_mask(i,j) = true;
                    boundary_val(i,j) = kbond*(psi(i,j)/xstp);
                end
            end
            % Upper boundary, inner points (i=1, 1<j<xnum)
            if(i==1 && j>1 && j<xn)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition: vx(i,j) = 0
                    % % equivalent to No slip vx=0: vx(i,j)-1/3*vx(i+1,j)=0.
                    Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                    Lf(invx,invx+2,-1/3*kbond); % Coefficient for vx(i+1,j)
                    %                     R(invx,1)          =   0;          % Right-hand-side part
                    boundary_mask(i,j) = true;
                elseif(bc==NEUMANN_BC || bc == PSI_BC)
                    % Neumann condition dvx/dy=0: vx(i+1,j)-vx(i,j)=0
                    %                 Same for psi boundary condition!
                    Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                    Lf(invx,invx+2,-1*kbond);    % Coefficient for vx(i+1,j)
                    %                     R(invx,1)           =   0;          % Right-hand-side part
                    boundary_mask(i,j) = true;
                end
            end
            % Lower boundary, inner points (i=ynum-1, 1<j<xnum)
            if(i==yn-1 && j>1 && j<xn)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet or No slip:
                    % % No slip vx=0: vx(i,j)-1/3*vx(i-1,j)=0
                    Lf(invx,invx,1*kbond); % Coefficient for vx(i,j)
                    Lf(invx,invx-2,-1/3*kbond); % Coefficient for vx(i-1,j)
                    %                     R(invx,1)         =   0; % Right part
                    boundary_mask(i,j) = true;
                elseif(bc==NEUMANN_BC || bc == PSI_BC)
                    % Neumann condition dvx/dy=0: vx(i,j)-vx(i-1,j)=0
                    %                 of psi condition
                    Lf(invx,invx,1*kbond); % Coefficient for vx(i,j)
                    Lf(invx,invx-2,-1*kbond); % Coefficient for vx(i-1,j)
                    %                     R(invx,1)           =   0; % Right-hand-side part
                    boundary_mask(i,j) = true;
                end
            end
            %Internal nodes: mu(d2vx/dx2+d2vx/dy2) + r d/dx(dvx/dx +
            %dvy/dy) = dP/dx + (r + mu + lambda)da/dx
        else
            %dSxx/dx=2*muc(i+1,j+1)*(vx(i,j+1)-vx(i,j))/dx^2-2*muc(i+1,j)*(vx(i,j)-vx(i,j-1))/dx^2
            Lf(invx,invx+yn*2,2*muc(i+1,j+1)/xstp^2);                         % Coefficient for vx(i,j+1)
            Lf(invx,invx-yn*2,2*muc(i+1,j)/xstp^2);                           % Coefficient for vx(i,j-1)
            Lf(invx,invx,-2*muc(i+1,j+1)/xstp^2-2*muc(i+1,j)/xstp^2);   % Coefficient for vx(i,j)
            
            %dSxy/dy=mug(i+1,j)*((vx(i+1,j)-vx(i,j))/dy^2+(vy(i+1,j)-vy(i+1,j-1))/dx/dy)-
            %         -mug(i,j)*((vx(i,j)-vx(i-1,j))/dy^2+(vy(i,j)-vy(i,j-1))/dx/dy)-
            Lf(invx, invx+2, mug(i+1,j)/ystp^2);                             % Coefficient for vx(i+1,j)
            Lf(invx, invx-2, mug(i,j)/ystp^2);                               % Coefficient for vx(i-1,j)
            Lf(invx, invx, -mug(i+1,j)/ystp^2-mug(i,j)/ystp^2); % ADD coefficient for vx(i,j)
            Lf(invx, invy+2, mug(i+1,j)/xstp/ystp);                          % Coefficient for vy(i+1,j)
            Lf(invx, invy+2-yn*2, -mug(i+1,j)/xstp/ystp);                         % Coefficient for vy(i+1,j-1)
            Lf(invx, invy, -mug(i,j)/xstp/ystp);                           % Coefficient for vy(i,j)
            Lf(invx, invy-yn*2, mug(i,j)/xstp/ystp);                            % Coefficient for vy(i,j-1)
            
            %d/dx(r*(dvx/dx + dvy/dy))
            Lf(invx, invx+yn*2, rc(i+1,j+1)/(xstp^2));            %Coeff. for vx(i,j+1)
            Lf(invx, invx, (-rc(i+1,j+1)-rc(i+1,j))/(xstp^2));           %Coeff. for vx(i,j)
            Lf(invx, invy+2, rc(i+1,j+1)/(xstp*ystp));             %Coerff. for vy(i+1,j)
            Lf(invx, invy, -rc(i+1,j+1)/(xstp*ystp));             %Coeff. for vy(i,j)
            Lf(invx, invx-yn*2, rc(i+1,j)/(xstp^2));            %Coeff. for vx(i,j-1)
            Lf(invx, invy+2-yn*2, -rc(i+1,j)/(xstp*ystp));      %Coeff. for vy(i+1,j-1)
            Lf(invx, invy-yn*2, rc(i+1,j)/(xstp*ystp));         %Coeff. for vy(i,j-1)
%             
%             Lf(invx, invx+yn*2, rx(i,j)/(xstp^2));            %Coeff. for vx(i,j+1)
%             Lf(invx, invx, (-rx(i,j)-rx(i,j))/(xstp^2));           %Coeff. for vx(i,j)
%             Lf(invx, invy+2, rx(i,j)/(xstp*ystp));             %Coerff. for vy(i+1,j)
%             Lf(invx, invy, -rx(i,j)/(xstp*ystp));             %Coeff. for vy(i,j)
%             Lf(invx, invx-yn*2, rx(i,j)/(xstp^2));            %Coeff. for vx(i,j-1)
%             Lf(invx, invy+2-yn*2, -rx(i,j)/(xstp*ystp));      %Coeff. for vy(i+1,j-1)
%             Lf(invx, invy-yn*2, rx(i,j)/(xstp*ystp));         %Coeff. for vy(i,j-1)
            
            
            %             Compute Right hand-side within the loop since it changes in
            %             every iteration depending on the value of p!!
            
            % -dP/dx=(P(i+1,j)-P(i+1,j+1))/dx
            %             Lf(invx, inp+3, kcont/xstp);                                     % Coefficient for P(i+1,j)
            %             Lf(invx, inp+3+yn*3, -kcont/xstp);                                    % Coefficient for P(i+1,j+1)
            % Right part:da/dx = (a(i+1,j+1)-a(i+1,j))/dx
            %             R(invx,1)               =   ((muc(i+1,j+1)+muc(i+1,j))/2 + (lambdac(i+1,j+1)+lambdac(i+1,j))/2) * (a(i+1,j+1)-a(i+1,j))/xstp;
        end
        
        % y-Stokes equation
        % Ghost vy unknowns (j=xnum) and boundary nodes (i=1, i=ynum, j=1, j=xnum-1)
        if(i==1 || i==yn || j==1 || j==xn-1 || j==xn)
            % Ghost vy unknowns (j=xnum: vy(i,j)=0
            if(j==xn)
                Lf(invy,invy,1*kbond);                    % Coefficient for vy(i,j)
                %                 R(invy,1)           =   0;
                boundary_mask(i,j) = true;
            end
            % Upper boundary, i=1
            if((i==1) && (j<xn))
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition vy = 0
                    Lf(invy,invy,1*kbond);
                    %                     R(invy,1) = 0;
                    boundary_mask(i,j) = true;
                elseif(bc==NEUMANN_BC)
                    % %            Neumann codition: -vy(i+1,j) + vy(i,j) = 0;
                    Lf(invy,invy,1*kbond);                    % Coefficient for vy(i,j)
                    Lf(invy,invy+2,-1*kbond);    %Coefficient for vy(i+1,j)
                    %                     R(invy,1)           =   0;
                    boundary_mask(i,j) = true;
                elseif(bc==PSI_BC)
                    %             psi boundary condtion:
                    Lf(invy,invy,1*kbond);
                    %                     R(invy,1) = kbond*(psi(i,j)/ystp);
                    boundary_mask(i,j) = true;
                    boundary_val = kbond*(psi(i,j)/ystp);
                end
            end
            %             Lower boundary, i=ynum
            if(i==yn && j < xn)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition vy = 0
                    Lf(invy,invy,1*kbond);
                    %                     R(invy,1) = 0;
                    boundary_mask(i,j) = true;
                elseif(bc==NEUMANN_BC)
                    % %           Neumann codition: vy(i,j) - vy(i-1,j) = 0;
                    Lf(invy,invy,1*kbond);   %Coeff. for vy(i,j)
                    Lf(invy, invy-2, -1*kbond); %Coeff. for vy(i-1,j)
                    %                     R(invy,1) = 0;               %Right side
                    boundary_mask(i,j) = true;
                elseif(bc==PSI_BC)
                    %                 psi boundary condtion:
                    Lf(invy,invy,1*kbond);
                    %                     R(invy,1) = kbond*(psi(i,j)/ystp);
                    boundary_mask(i,j) = true;
                    boundary_val(i,j) = kbond*(psi(i,j)/ystp);
                end
            end
            
            % Left boundary, inner points (j=1, 1<i<ynum)
            if(j==1 && i>1 && i<yn)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet or eqv. to no slip:
                    %             % No slip vy=0: vy(i,j)-1/3*vy(i,j+1)=0
                    Lf(invy,invy,1*kbond); % Coefficient for vy(i,j)
                    Lf(invy,invy+yn*2,-1/3*kbond); % Coefficient for vy(i,j+1)
                    %                     R(invy,1)=0;
                    boundary_mask(i,j) = true;
                elseif(bc==NEUMANN_BC || bc == PSI_BC)
                    % Neumann dvy/dx=0: vy(i,j)-vy(i,j+1)=0
                    %                 eqv. for psi boundary condtion
                    Lf(invy,invy,1*kbond);                   % Coefficient for vy(i,j)
                    Lf(invy,invy+yn*2,-1*kbond);                   % Coefficient for vy(i,j+1)
                    %                     R(invy,1)=0;
                    boundary_mask(i,j) = true;
                end
            end
            % Right boundary, inner points (j=xnum-1, 1<i<ynum)
            if(j==xn-1 && i>1 && i<yn)
                if(bc==DIRICHLET_BC)
                    %                 Dirichlet or no slip: vy=0: vy(i,j)-1/3*vy(i,j-1)=0
                    Lf(invy,invy,1*kbond); % Coefficient for vy(i,j)
                    Lf(invy,invy-yn*2,-1/3*kbond); % Coefficient for vy(i,j-1)
                    %                     R(invy,1)=0;
                    boundary_mask(i,j) = true;
                elseif(bc==NEUMANN_BC || bc==PSI_BC)
                    % Free slip dvy/dx=0: vy(i,j)-vy(i,j-1)=0
                    %                 eqv. to psi boundary condition
                    Lf(invy,invy,1*kbond);                    % Coefficient for vy(i,j)
                    Lf(invy,invy-yn*2,-1*kbond);                   % Coefficient for vy(i,j-1)
                    %                     R(invy,1)=0;
                    boundary_mask(i,j) = true;
                end
            end
            %Internal nodes: mu(d2vy/dx2+d2vy/dy2)-dP/dy=-gy*RHO
        else
            %dSyy/dy=2*muc(i+1,j+1)*(vy(i+1,j)-vy(i,j))/dy^2-2*muc(i,j+1)*(vy(i,j)-vy(i-1,j))/dy^2
            Lf(invy, invy+2, 2*muc(i+1,j+1)/ystp^2);                 % Coefficient for vy(i+1,j)
            Lf(invy, invy-2, 2*muc(i,j+1)/ystp^2);                   % Coefficient for vy(i-1,j)
            Lf(invy, invy, -2*muc(i+1,j+1)/ystp^2-2*muc(i,j+1)/ystp^2); % Coefficient for vy(i,j)
            
            %dSxy/dx=mus(i,j+1)*((vy(i,j+1)-vy(i,j))/dx^2+(vx(i,j+1)-vx(i-1,j+1))/dx/dy)-
            %         -mus(i,j)*((vy(i,j)-vy(i,j-1))/dx^2+(vx(i,j)-vx(i-1,j))/dx/dy)-
            Lf(invy,invy+yn*2, mug(i,j+1)/xstp^2);                     % Coefficient for vy(i,j+1)
            Lf(invy,invy-yn*2, mug(i,j)/xstp^2);                       % Coefficient for vy(i,j-1)
            Lf(invy,invy, -mug(i,j+1)/xstp^2-mug(i,j)/xstp^2); % ADD coefficient for vy(i,j)
            Lf(invy,invx+yn*2, mug(i,j+1)/xstp/ystp);                  % Coefficient for vx(i,j+1)
            Lf(invy,invx+yn*2-2, -mug(i,j+1)/xstp/ystp);                 % Coefficient for vx(i-1,j+1)
            Lf(invy,invx, -mug(i,j)/xstp/ystp);                   % Coefficient for vx(i,j)
            Lf(invy,invx-2, mug(i,j)/xstp/ystp);                    % Coefficient for vx(i-1,j)
            
            
            %d/dy(r*(dvx/dx + dvy/dy))
            Lf(invx, invx+yn*2, rc(i+1,j+1)/(xstp*ystp));         %Coeff. for vx(i,j+1)
            Lf(invx, invx, -rc(i+1,j+1)/(xstp*ystp));             %Coeff. for vx(i,j)
            Lf(invx, invy+2, rc(i+1,j+1)/(ystp^2));               %Coerff. for vy(i+1,j)
            Lf(invx, invy, (-rc(i+1,j+1) - rc(i,j+1))/(ystp^2));              %Coeff. for vy(i,j)
            Lf(invx, invx-2+yn*2, -rc(i,j+1)/(xstp*ystp));      %Coeff. for vx(i-1,j+1)
            Lf(invx, invx-2, rc(i,j+1)/(xstp*ystp));            %Coeff. for vx(i-1,j)
            Lf(invx, invy-2, rc(i,j+1)/(ystp^2));               %Coeff. for vy(i-1,j)
            
            
%             Lf(invx, invx+yn*2, ry(i,j)/(xstp*ystp));         %Coeff. for vx(i,j+1)
%             Lf(invx, invx, -ry(i,j)/(xstp*ystp));             %Coeff. for vx(i,j)
%             Lf(invx, invy+2, ry(i,j)/(ystp^2));               %Coerff. for vy(i+1,j)
%             Lf(invx, invy, (-ry(i,j) - ry(i,j))/(ystp^2));              %Coeff. for vy(i,j)
%             Lf(invx, invx-2+yn*2, -ry(i,j)/(xstp*ystp));      %Coeff. for vx(i-1,j+1)
%             Lf(invx, invx-2, ry(i,j)/(xstp*ystp));            %Coeff. for vx(i-1,j)
%             Lf(invx, invy-2, ry(i,j)/(ystp^2));               %Coeff. for vy(i-1,j)
            
            %             RHS needs modification in every step depending on pressure
            %             values!!
            % -dP/dy=(P(i,j+1)-P(i+1,j+1))/dx
            %             Lf(invy,inp+yn*3, kcont/ystp);                             % Coefficient for P(i,j+1)
            %             Lf(invy,inp+3+yn*3, -kcont/ystp);                            % Coefficient for P(i+1,j+1)
            %             % Right part: da/dy = (a(i+1,j+1)-a(i,j+1))/dy
            %             R(invy,1)               = ((muc(i+1,j+1)+muc(i,j+1))/2 + (lambdac(i+1,j+1)+lambdac(i,j+1))/2)*(a(i+1,j+1)-a(i,j+1))/ystp;
        end
        
    end
end



% % Constraints when neumann boundary for displacement is used for all the boundaries
% if (bc==NEUMANN_BC)
%     % For rigid body translation
%     % % Let's set the average of the each of the displacement components to be
%     % % zero!
%     for i = 2:yn-1
%         for j = 2:xn-1
%             invx = ((j-1)*yn+i)*2 - 1;
%             invy = invx + 1;
%             Lf(1+xn*yn*2,invx, 1*kbond); % Coefficient for vx(i+1,j)
%             Lf(2+xn*yn*2,invy, 1*kbond); % Coefficient for vx(i-1,j) %NOT
% %             probably a mistake, index in comment does not match the one in the statement!!
%         end
%     end
%     R(1+xn*yn*2,1) = 0;
%     R(2+xn*yn*2,1) = 0;
%
%
%     % For rigid body rotation
%     % At all interior points, mean of rigid body rotation tensor components set
%     % to zero:
%     % i.e. dvx/dy - dvy/dx = 0.
%     % dvx/dy: (vx(i+1,j)-vx(i-1,j))/(2*ystp) - (vy(i,j+1)-vy(i,j-1))/(2*ystp) = 0;
%     for i = 2:yn-1
%         for j = 2:xn-1
%             invx = ((j-1)*yn+i)*2 - 1;
%             invy = invx + 1;
%             Lf(3+xn*yn*2,invx + 2, kbond/(2*xstp)); % Coefficient for vx(i+1,j)
%             Lf(3+xn*yn*2,invx - 2, -kbond/(2*xstp)); % Coefficient for vx(i-1,j)
%
%             Lf(3+xn*yn*2,invy + yn*2, -kbond/ystp); % Coefficient for vy(i,j+1)
%             Lf(3+xn*yn*2,invy - yn*2, kbond/ystp); % Coefficient for vy(i,j-1)
%         end
%     end
%     R(3+xn*yn*2,1) = 0;
%     ecnstrnt_indx = 4;
%
%     % If single point is considered instead of for all interior points:
%     % Fix the center, i.e. the vx and vy corresponding to (ynum/2,xnum/2):
%     % i.e. vx(ynum/2,xnum/2) = 0;
%     % and vy(ynum/2,xnum/2) = 0;
%
%     % top left corner
%     % i = 2;
%     % j = 2;
%
%     % top middle
%     % i = 2;
%     % j = ceil(xnum/2);
%
%     % top right corner
%     % i = 2;
%     % j = xnum-1;
%
%     % center
%     % i = ceil(ynum/2);
%     % j = ceil(xnum/2);
%
%     % bottom left corner
%     % i = ynum - 1;
%     % j = 2;
%
%     % bottom middle
%     % i = ynum - 1;
%     % j = ceil(xnum/2);
%
%     % bottom right
%     % i = ynum - 1;
%     % j = xnum - 1;
%
%     % invx = ((j-1)*ynum+i)*3-1;
%     % invy = invx + 1;
%     %
%     % % % vx(i,j) = 0;
%     % L(1+xnum*ynum*3,invx) = 1*kbond;
%     % R(1+xnum*ynum*3,1) = 0;
%     % %
%     % % % vy(i,j) = 0;
%     % L(2+xnum*ynum*3,invy) = 1*kbond;
%     % R(2+xnum*ynum*3,1) = 0;
%
%
%
%     % Now set dvx/dy - dvy/dx = 0 as another constraint!
%     % dvx/dy: (vx(i+1,j)-vx(i-1,j))/(2*ystp) - (vy(i,j+1)-vy(i,j-1))/(2*ystp) = 0;
%     % At single point
%     % L(3+xnum*ynum*3,invx + 3       )    =    kcont/(2*xstp); % Coefficient for vx(i+1,j)
%     % L(3+xnum*ynum*3,invx - 3)    =   -kcont/(2*xstp); % Coefficient for vx(i-1,j)
%     % % % % %
%     % % % % %-dvy/dx: -(vy(i,j+1)+vy(i,j-1))/(2*xstp)
%     % L(3+xnum*ynum*3,invy + ynum*3  )    =   -kcont/ystp; % Coefficient for vy(i,j+1)
%     % L(3+xnum*ynum*3,invy - ynum*3)    =   kcont/ystp; % Coefficient for vy(i,j-1)
%     %
%     % R(3+xnum*ynum*3,1) = 0;
%
% %     if (pChoice == ENFORCE_CC) %Since CC is enforced, extra constraint needed here to fix pressure value to make it full rank!
% %         inp     =   ((p0cellx-1)*yn+p0celly)*3-2;
% %         Lf(ecnstrnt_indx + xn*yn*3,inp, 1*kbond);  %Coeff. for p0cell
% %         R(ecnstrnt_indx + xn*yn*3,1) = p0cell;
% %         ecnstrnt_indx = ecnstrnt_indx + 1;
% %     end
% end

% if (bc ~= NEUMANN_BC && pChoice == ENFORCE_CC)
%     inp     =   ((p0cellx-1)*yn+p0celly)*3-2;
%     Lf(1+xn*yn*3,inp, 1*kbond);  %Coeff. for p0cell
%     R(1+xn*yn*3,1) = p0cell;
%     ecnstrnt_indx = ecnstrnt_indx + 1;
% end


% Create sparse matrix from the indices and values
L = sparse(i_indx(1:vals_cntr),j_indx(1:vals_cntr),L_vals(1:vals_cntr),ecnstrnt+xn*yn*2,xn*yn*2);

toc;
display('matrix built');
% save(['staggeredAug2D_' num2str(xn) '_.mat'],'L','xn','yn'); 
else
    load staggeredAug2D_2000_.mat
    display('matrix read');
end

% ahilelai talako purai comment:
% % write the matrix to solve with PetsC:
% writeAandBtoFile(L,R,'petsC/A.txt','petsC/b.txt');
% tic;
%
% Augmented Lagrangian Iteration loop:
% Initialize p:
p = zeros(yn,xn);
v = zeros(yn,xn,2);
max_it = 200;
fl = zeros(max_it,1);
%     tol = 1e-12; maxit_solve = 40;
err = zeros(max_it,1);
momentum_residue = zeros(max_it,1);
for indx = 1:max_it
    rhs = setRhs(p);
%     S = L\rhs;
    [LL,UU] = ilu(L,struct('type','ilutp','droptol',1e-5));
%     [S,fl(indx),rr1,it1,rv1] = bicgstab(L,rhs,tol,maxit_solve,LL,UU);
    [S,fl(indx),rr1,it1,rv1] = bicgstab(L,rhs,1e-5,20,LL,UU);
% [S fl(indx)] = bicgstab(L,rhs,1e-5,100);
momentum_residue(indx) = norm(L*S-rhs);
% [S fl1] = gmres(L,rhs);
%     Rearrange the solution in matrix form
    for i=1:1:yn
        for j=1:1:xn
        % Global index for vx, vy in S()
            invx=((j-1)*yn+i)*2-1; % P
            invy=invx+1;
            % vx
            v(i,j,1)=S(invx);
            % vy
            v(i,j,2)=S(invy);
        end
    end
    div_v = div2D_new(v,xstp,ystp);
    curr_err = div_v(1:end-1,1:end-1)+a(2:end,2:end);
    p(2:end,2:end) = p(2:end,2:end) - rc(2:end,2:end).*curr_err;
    err(indx) = norm(curr_err);
%     p(2:end,2:end) = p(2:end,2:end) - r*div_v(1:end-1,1:end-1) - r*a(2:end,2:end);
end
stem(err);
figure, stem(momentum_residue), title('momentum residue');
figure,
subplot(121), imagesc(div_v(1:end-1,1:end-1)), axis image, subplot(122), imagesc(a(2:end,2:end)), axis image;
figure,
imagesc(div_v(1:end-1,1:end-1)+a(2:end,2:end)), title('constraint');
vx = v(:,:,1);
vy = v(:,:,2);
vx1=zeros(yn,xn);
vy1=zeros(yn,xn);
% Process internal Grid points
for i=2:1:yn-1
    for j=2:1:xn-1
        % vx
        vx1(i,j)=(vx(i-1,j)+vx(i,j))/2;
        % vy
        vy1(i,j)=(vy(i,j-1)+vy(i,j))/2;
    end
end

figure,
imagesc(p), hold on, quiver(vx1,vy1), axis ij image;
kk= 1;
    
% %Obtaining vector of solutions S()
% % S=L\R;
%
% tol = 1e-12; maxit = 40;
% [LL,UU] = ilu(L,struct('type','ilutp','droptol',1e-5));
% [S,fl1,rr1,it1,rv1] = bicgstab(L,R,tol,maxit,LL,UU);
%
%
% % options = AMGinit(L);
% % [prec options] = AMGfactor(L,options);
% % [S, options] = AMGsolver(S,prec,options,R);
%
%
% toc;
% display('system solver timed');
% % Reload solutions to 2D p(), vx(), vy() arrays
% % Dimensions of arrays are reduced compared to the basic grid
% p=zeros(yn,xn);
% vy=zeros(yn,xn);
% vx=zeros(yn,xn);
% % Process all Grid points
% for i=1:1:yn
%     for j=1:1:xn
%         % Global index for P, vx, vy in S()
%         inp=((j-1)*yn+i)*3-2; % P
%         invx=inp+1;
%         invy=inp+2;
%         % P
%         p(i,j)=S(inp)*kcont;
%         % vx
%         vx(i,j)=S(invx);
%         % vy
%         vy(i,j)=S(invy);
%     end
% end
%
% save('results2D.mat','vx','vy','p','L','R');
%
% % Compute vx,vy for internal nodes
% % Note here too, vx1 and vy1 will have dummy 1st row and 1st col with zero!
% % So compatible with the dummy atrophy!
% vx1=zeros(yn,xn);
% vy1=zeros(yn,xn);
% % Process internal Grid points
% for i=2:1:yn-1
%     for j=2:1:xn-1
%         % vx
%         vx1(i,j)=(vx(i-1,j)+vx(i,j))/2;
%         % vy
%         vy1(i,j)=(vy(i,j-1)+vy(i,j))/2;
%     end
% end
%
%
%
% %Plotting solution
% % Making new figure
% figure;
%
% % Plotting pressure as colormap
% % subplot(1,2,1);
% % pcolor(xc/1000,yc/1000,p(2:1:ynum,2:1:xnum));%*1e-9);      % making a colormap
% % shading interp;     % making smooth transitions between colors
% imagesc(p(2:1:yn,2:1:xn));%*1e-9);
% colorbar;           % showing a colorbar for the map
% hold on;            % continuing plotting on the colormap
% % Plotting velocity vector as arrows using internal nodes only
% % quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
% % quiver(vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
% quiver(vx1(2:1:yn-1,2:1:xn-1)*10000,vy1(2:1:yn-1,2:1:xn-1)*10000,'k','Autoscale','Off'); % making field of arrows
% hold off;           % stop plotting on the colormap
% box on;             % making a box around the plot
% title('Pressure (color), velocity (arrows)'); % title for the plot
% % xlabel('x, m');        % title for the horizontal axis
% % ylabel('y, m');        % title for the vertical axis
% axis ij image ;     % directing vertical axis downward, making proper dimensions
% % axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
%
%
% % Plotting viscosity
% figure,
% imagesc(mug), title(['muBrain: ' num2str(muBrain) ' muCSF: ' num2str(muCSF)]);
% axis equal;
% % Plotting density as colormap
% % subplot(1,2,2);
% % pcolor(x/1000,y/1000,rho);      % making a colormap
% % shading interp;     % making smooth transitions between colors
% % colorbar;           % showing a colorbar for the map
% % hold on;            % continuing plotting on the colormap
% % Plotting velocity vector as arrows using internal nodes only
% % quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
% % hold off;           % stop plotting on the colormap
% % box on;             % making a box around the plot
% % title('Density (color, kg/m^3), velocity (arrows)');   % title for the plot
% % xlabel('x, km');        % title for the horizontal axis
% % ylabel('y, km');        % title for the vertical axis
% % axis ij image ;     % directing vertical axis downward, making proper dimensions
% % axis([0 xsize/1000 0 ysize/1000]); % Making axes limits
%
% figure,
% subplot(121),
% % Display atrophy except the dummy ones (i.e. 1st row and 1st column)
% imagesc(a(2:end,2:end)), title('atrophy, a');
% % axis image;
% axis equal;
%
% soltn(:,:,1) = vx;
% soltn(:,:,2) = vy;
% div_ans = div2D_new(soltn,xstp,ystp);
% % Last row/col contains false div. so remove it
% div_ans = div_ans(1:end-1,1:end-1);
% % Now add dummy on first row and col to make it comparible with a:
% div_ans = [zeros(1,xn); zeros(yn-1,1) div_ans];
% % Display divergence of the solution except for the dummy row and col.
% subplot(122), imagesc(div_ans(2:end,2:end)), title('divergence of the solution displ. field ');
% axis equal;
%
%
% % clims = [min([min(a) min(div_ans)]) max([max(a) max(div_ans)])];
% figure,
% % subplot(121), imagesc(a,clims), title('atrophy');
% % subplot(122), imagesc(-1*(div_ans),clims), title('divergence');
% % subplot(121), imagesc(a), title('atrophy');
% % subplot(122), imagesc(div_ans), title('divergence');
% imagesc(a+div_ans), title('constraint');
% axis equal;
%
% % Re-add the removed part with zero velocities.
% % Note that in using vx1, vx1 already has a zero first row and col due to
% % dummy addition during interpolations from vx before this concatenation.
% if(USE_REAL_DATA == true)
%     [m n] = size(div_slice);
%     vc = size(vx1,2);
%     vx2 = cat(1,zeros(min_max_row_col(1)-extra_CSF-1-1, vc),vx1,zeros(m-min_max_row_col(2)-extra_CSF-1+1,vc));
%     vx2 = cat(2,zeros(m, min_max_row_col(3)-extra_CSF-1-1),vx2,zeros(m, n-min_max_row_col(4)-extra_CSF-1+1));
%
%     vy2 = cat(1,zeros(min_max_row_col(1)-extra_CSF-1-1, vc),vy1,zeros(m-min_max_row_col(2)-extra_CSF-1+1,vc));
%     vy2 = cat(2,zeros(m, min_max_row_col(3)-extra_CSF-1-1),vy2,zeros(m, n-min_max_row_col(4)-extra_CSF-1+1));
%
%     v2(:,:,1) = vx2;
%     v2(:,:,2) = vy2;
% end

end