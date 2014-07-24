% 3D Solution of the Brain deformation with irregular domain (or brain-CSF
% considered as an interface and extending the boundary to the skull.
% Regular staggered grid using a pressure-velocity formulation
% for a medium with varying lame's parameters.

function AdLemTaras3D(seg_mask,div_map,USE_REAL_DATA,muBrain,muRatio,lambdaBrain,lambdaRatio,voi)

if(nargin < 3)
    USE_REAL_DATA = false;
elseif (nargin < 7)
    if(USE_REAL_DATA == true);
        error('too few arguments for real data case');
    end
end

% Since mu and lambda Brain are greater than CSF, the ratio muB/muC given
% as input in the form muRatio should be greater or equal to 1.
if(lambdaRatio < 1 || muRatio < 1)
    error('invalid muRatio or invalid lamdaRatio, they should be greater than one');
end

if (USE_REAL_DATA == true)
    %     [brain_mask min_max_row_col] = getRectConvexHull(seg_mask);
    % or, extract a small 3D volume
    top_left = voi(1,:);  %  %top-left
    vol_size = voi(2,:);  %     %vol-size
    brain_mask  = seg_mask(top_left(1):top_left(1)+vol_size(1),...
        top_left(2):top_left(2)+vol_size(2),top_left(3):top_left(3)+vol_size(3));
    
    a = div_map(top_left(1):top_left(1)+vol_size(1),...
        top_left(2):top_left(2)+vol_size(2),top_left(3):top_left(3)+vol_size(3));
    
    [yn xn zn] = size(brain_mask);
    %pad with extra_CSF width of zeros all around.
    extra_CSF = 2;
    brain_mask = padarray(brain_mask,[extra_CSF extra_CSF extra_CSF]);
    
    % distributed source in the padded area excluding
    % the 12 edges where CC is not enforced! i.e. ((xn+2w)(yn+2w)(zn+2w) -
    % xn.yn.zn - (4(xn-3 + yn-3 + zn-3)))
    pix_considered = (xn+2*extra_CSF)*(yn+2*extra_CSF)*(zn+2*extra_CSF) - ...
        xn*yn*zn - (4*(xn-3 + yn-3 + zn-3));
    
    a = padarray(a,[extra_CSF extra_CSF extra_CSF], -sum(a(:))/pix_considered);
    
    %Grid nodes have one dimension greater than the cell centers!
    yn = yn+extra_CSF*2+1;
    xn = xn+extra_CSF*2+1;
    zn = zn+extra_CSF*2+1;
    
    % add a dummy atrophy corresponding to dummy pressure values!
    a = cat(3,zeros(yn-1,xn-1,1),a);
    a = [zeros(1,xn,zn); zeros(yn-1,1,zn) a];
    
    % Model size from voxel resolution : xres X yres X zres
    % For 1mm X 1mm X 1mm
    xres = 1e-3; yres = 1e-3; zres = 1e-3;
    xsize   =   (xn-1)*xres;        % Horizontal
    ysize   =   (yn-1)*yres;        % Vertical
    zsize   =   (zn-1)*zres;        % Anterior-Posterior
    
else
    
    
    % Numerical model parameters
    % Model size, m
    xsize   =   1;        % Horizontal rightwards
    ysize   =   1;        % Vertical
    zsize   =   1;        % Horizontal towards the screen
    
    % Numbers of nodes
    xn    =   5;             % Horizontal
    yn    =   5;             % Vertical
    zn    =   5;             % Horizontal towards the screen
    
    % Call a function:
    center_mask = [round(xn/2) round(yn/2) round(zn/2)];
    % brain_mask = sphereMask(center_mask,xnum/4,xnum-1,ynum-1);
    % brain_mask = ellipsoidMask(center_mask,xnum/3,ynum/4,znum/4,xnum-1,ynum-1,znum-1);
    % brain_mask = ellipsoidMask(center_mask,7,4,4,xnum-1,ynum-1,znum-1);
    brain_mask = sphereMask(center_mask,floor(xn/4),xn-1,yn-1,zn-1);
    % brain_mask = sphereMask(center_mask,7,xnum-1,ynum-1,znum-1);
    
    % Get the atrophy data at the cell-centers
    % atrophy at the center of the cells
    a = computeAtrophy3D(yn-1,xn-1,zn-1);
    
    % distributed source outside brain, excluding 12 edges where CC is not
    % enforced!
    a(brain_mask==false) = -sum(a(:))/(xn*yn*zn - (4*(xn-3 + yn-3 + zn-3)) - sum(brain_mask(:)));
    % add a dummy atrophy corresponding to dummy pressure values!
    a = cat(3,zeros(yn-1,xn-1,1),a);
    a = [zeros(1,xn,zn); zeros(yn-1,1,zn) a];
    
end

NS     =  xn*yn*zn*4;     % Total number of variables to compute values of.

% Grid step
xstp    =   xsize/(xn-1); % Horizontal
ystp    =   ysize/(yn-1); % Vertical
zstp    =   zsize/(zn-1); % Away from the screen
% % Create vectors for nodal points positions (basic nodes)
% x       =   0:xstp:xsize;   % Horizontal
% y       =   0:ystp:ysize;   % Vertical
% z       =   0:zstp:zsize;   %Horizontal towards the screen

% Create vectors for cell centers positions (staggered nodes)
% xc      =   xstp/2:xstp:xsize-xstp/2; % Horizontal
% yc      =   ystp/2:ystp:ysize-ystp/2; % Vertical
% zc      =   zstp/2:ystp:zsize-zstp/2; % Horizontal towards the

% Model viscosity

%     muBrain = 2500;
%     lambdaBrain = 2500;

muCSF = muBrain/muRatio;
lambdaCSF = lambdaBrain/lambdaRatio;


% mu and lambda in center of the cells!
muc = zeros(yn-1,xn-1,zn-1);
lambdac = zeros(yn-1,xn-1,zn-1);

muc(brain_mask==true) = muBrain;
muc(brain_mask==false) = muCSF;
lambdac(brain_mask==true) = lambdaBrain;
lambdac(brain_mask==false) = lambdaCSF;

% Interpolate mu at different edges from values at centers of the cells
muxy = interpolateAtEdgeFromCellCenter(muc,3);
muxz = interpolateAtEdgeFromCellCenter(muc,1);
muyz = interpolateAtEdgeFromCellCenter(muc,2);

% PACHHI HATAAUNU PARCHHA YO, EUTAI NAAM RAAKHNE!
muyx = muxy; muzx = muxz;

% Now Add dummy zeros on 1st row,col and plane of centered cells:
muc = cat(3,zeros(yn-1,xn-1,1),muc);
muc = [zeros(1,xn,zn); zeros(yn-1,1,zn) muc];

% Interpolate mu at different edges from values at centers of the cells
% lambdaxy = interpolateAtEdgeFromCellCenter(lambdac,3);
% lambdaxz = interpolateAtEdgeFromCellCenter(lambdac,1);
% lambdayz = interpolateAtEdgeFromCellCenter(lambdac,2);

% Now Add dummy zeros on 1st row,col and plane of centered cells:
lambdac = cat(3,zeros(yn-1,xn-1,1),lambdac);
lambdac = [zeros(1,xn,zn); zeros(yn-1,1,zn) lambdac];

% Computing Kcont and Kbond coefficients
kcont   =   1; %2*muBrain/(xstp+ystp+zstp);
kbond   =   1; %4*muBrain/(xstp+ystp+zstp)^2;
% kcont   =   2*muCSF/(xstp+ystp+zstp);
% kbond   =   4*muCSF/(xstp+ystp+zstp)^2;

DIRICHLET_BC = 0;
NEUMANN_BC = 1;

RELAX_CC = 1;
ENFORCE_CC = 0;

% Boundary Condition Options:
bc = DIRICHLET_BC; %dirichlet condition on displacement.
% bc = NEUMANN_BC; %neumann condition on displacement.

% Pressure condition in one cell (i,j)
p0cell  =   20;

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

% Number of Extra constraints required, 3 for translation, 3 for rotation
% if (bc == NEUMANN_BC)
if (bc == NEUMANN_BC)
    ecnstrnt = 6;
else
    ecnstrnt = 0;
end

if(pChoice == ENFORCE_CC)   %if CC enforced, extra constraints on pressure required to make full rank matrix.
    ecnstrnt = ecnstrnt + 1;
end

% index keeping track of the ecnstrnt
ecnstrnt_indx = 1;


% Top left:
p0cellx = 3; %indx j
p0celly = 3; %indx i
p0cellz = 2; %indx k

% Top right:
% p0cellx = xnum-1; %indx j
% p0celly = 2; %indx i
% p0cellz = 2; %indx k

% Bottom left:
% p0cellx = 3; %indx j
% p0celly = ynum-1; %indx i
% p0cellz = 2; %indx k

% Center
% p0cellx = round(xnum/2); %indx j
% p0celly = round(ynum/2); %indx i
% p0cellz = round(znum/2); %indx k

% Bottom right
% p0cellx = xn-2; %indx j
% p0celly = yn-2; %indx i
% p0cellz = zn-1; %indx k

% Arbitrary point:
% p0cellx = round(xnum/2)+4;
% p0celly = round(ynum/2)-4;
% p0cellz = round(znum/2)-4;


% Matrix of coefficients initialization
% L       =   sparse(ecnstrnt+ NS, NS);
% L       =   zeros(ecnstrnt+ NS, NS);

% Tentative number of non-zeros in L matrix!
% NEED TO BE CALCULATED AND UPDATED for 3D
tot_non_zeros = 4*(xn+yn) + 4*(xn-1)*(yn-1) + 12*((yn-3)*(xn-2) + ...
    (xn-3)*(yn-2)) + 5 + 2*yn + 2*xn; %last ko 2*yn+2*xn ettikai extra
% global i_indx;
i_indx = zeros(tot_non_zeros,1);
% global j_indx;
j_indx = zeros(tot_non_zeros,1);
% global L_vals ;
L_vals = zeros(tot_non_zeros,1);


% Vector of right part initialization
R       =   zeros(ecnstrnt+ NS,1);

% global vals_cntr;
vals_cntr = 0;

    function Lf(ii,jj,xx)
        %     global vals_cntr;
        %     global i_indx;
        %     global j_indx;
        %     global L_vals;
        
        vals_cntr = vals_cntr + 1;
        i_indx(vals_cntr) = ii;
        j_indx(vals_cntr) = jj;
        L_vals(vals_cntr) = xx;
    end


% Solving x-Stokes, y-Stokes and continuity equations
% x-Stokes: mu(d2vx/dx2+d2vx/dy2)-dP/dx=(mu+lambda)da/dx
% y-Stokes: mu(d2vy/dx2+d2vy/dy2)-dP/dy=gy*RHO + (mu+lambda)da/dy
% continuity: dvx/dx+dvy/dy=-a
% Compose matrix of coefficients L()
% and vector (column) of right parts R()
% Boundary conditions: free slip
% Process all Grid points

tic;
for i=1:1:yn
    for j=1:1:xn
        for k = 1:1:zn
            
            % Global index for P, vx, vy and vz in the current node
            inp     =   ((k-1)*xn*yn + (j-1)*yn+i)*4-3; % P
            invx    =   inp+1;
            invy    =   inp+2;
            invz    =   inp+3;


            %         It's the same for atrophy too, so we can use this index to access
            %         atrophy values too!
            
            % Continuity equation
            % Ghost pressure unknowns (i=1, j=1, k = 1) and boundary nodes (12
            % edges + one cell)
            
            if ( (i==1) || (j==1) || (k==1) || (i==2 && j==2 && k>1) || (i==yn && j==2 && k>1) || (i==2 && j==xn && k>1) || (i==yn && j==xn && k>1) ...
                    || (i==2 && k==2 && (2<j<xn)) || (i==yn && k==2 && (2<j<xn)) || (i==2 && k==zn && (2<j<xn)) || (i==yn && k==zn && (2<j<xn)) ...
                    || ((2<i<yn) && j==2 && k==2) || ((2<i<yn) && j==2 && k==zn) || ((2<i<yn) && j==xn && k==2) || ((2<i<yn) && j==xn && k==zn) ...
                    || (i==p0celly && j==p0cellx && k==p0cellz && pChoice == RELAX_CC))
                % Ghost pressure unknowns (i=1, j=1, k=1): P(i,j)=0
                if(i==1 || j==1 || k==1)
                    Lf(inp,inp, 1*kbond);    % Coefficient for P(i,j,k)
                    R(inp,1  )          =   0;          % Right-hand-side part
                end
                
                %Upper and lower horizontal left edges dP/dx=0 => P(i,j,k)-P(i,j+1,k)=0
                if((i==2 && j==2 && k>1) || (i==yn && j==2 && k>1))
                    Lf(inp,inp, kbond);   % Coefficient for P(i,j,k)
                    Lf(inp,inp+yn*4, -1*kbond);   % Coefficient for P(i,j+1,k)
                    R(inp,1)            =   0;          % Right-hand-side part
                end
                
                % Upper and lower horizontal right edges dP/dx=0 => P(i,j,k)-P(i,j-1,k)=0
                if((i==2 && j==xn && k>1) || (i==yn && j==xn && k>1))
                    Lf(inp,inp, 1*kbond);   % Coefficient for P(i,j,k)
                    Lf(inp,inp-yn*4, -1*kbond);   % Coefficient for P(i,j-1,k)
                    R(inp,1)            =   0;          % Right-hand-side part
                end
                
                % Upper and lower horizontal front edges (internal points)
                %             dP/dz=0 => P(i,j,k)-P(i,j,k+1)=0
                if((i==2 && (2<j<xn) && k==2) || (i==yn && (2<j<xn) && k==2))
                    Lf(inp,inp, 1*kbond);   % Coefficient for P(i,j,k)
                    Lf(inp,inp+xn*yn*4, -1*kbond);   % Coefficient for P(i,j,k+1)
                    R(inp,1)            =   0;          % Right-hand-side part
                end
                % Upper and lower horizontal back edges dP/dz=0 =>
                % P(i,j,k)-P(i,j,k-1)=0
                if((i==2 && (2<j<xn) && k==zn) || (i==yn && (2<j<xn) && k==zn))
                    Lf(inp,inp, 1*kbond);   % Coefficient for P(i,j,k)
                    Lf(inp,inp-xn*yn*4, -1*kbond);   % Coefficient for P(i,j,k-1)
                    R(inp,1)            =   0;          % Right-hand-side part
                end
                
                % Front and Back vertical left edges (internal points)
                %             dP/dx=0 => P(i,j,k)-P(i,j+1,k)=0
                if(((2<i<yn) && j==2 && k==2) || ((2<i<yn) && j==2 && k==zn))
                    Lf(inp,inp, 1*kbond);   % Coefficient for P(i,j,k)
                    Lf(inp,inp+yn*4, -1*kbond);   % Coefficient for P(i,j+1,k)
                    R(inp,1)            =   0;          % Right-hand-side part
                end
                % Front and Back vertical right edges dP/dx=0 => P(i,j,k)-P(i,j-1,k)=0
                if(((2<i<yn) && j==xn && k==2) || ((2<i<yn) && j==xn && k==zn))
                    Lf(inp,inp, 1*kbond);   % Coefficient for P(i,j,k)
                    Lf(inp,inp-yn*4, -1*kbond);   % Coefficient for P(i,j-1,k)
                    R(inp,1)            =   0;          % Right-hand-side part
                end
                
                % One cell
                if (i==p0celly && j==p0cellx && k==p0cellz && pChoice == RELAX_CC) %relax CC at p0cell
                    Lf(inp,inp, 1*kbond);    % Coefficient for P(i,j)
                    R(inp,1)            =   p0cell;     % Right-hand-side part
                end
                %
                %Internal nodes: dvx/dx + dvy/dy + dvz/dz = -a
            else %corresponds to divergence at: p(i,j,k) & a(i,j,k)
                %dvx/dx=(vx(i-1,j,k-1)-vx(i-1,j-1,k-1))/dx
                Lf(inp,invx - 4 - xn*yn*4, kcont/xstp); % Coefficient for vx(i-1,j,k-1)
                Lf(inp,invx - 4 - yn*4 - xn*yn*4, -kcont/xstp); % Coefficient for vx(i-1,j-1,k-1)
                
                %dvy/dy=(vy(i,j-1,k-1)-vy(i-1,j-1,k-1))/dy
                Lf(inp,invy - yn*4 - xn*yn*4, kcont/ystp); % Coefficient for vy(i,j-1,k-1)
                Lf(inp,invy - 4 - yn*4 - xn*yn*4, -kcont/ystp); % Coefficient for vy(i-1,j-1,k-1)
                
                %dvz/dz=(vz(i-1,j-1,k)-vz(i-1,j-1,k-1))/dz
                Lf(inp,invz - 4 - yn*4, kcont/zstp); % Coefficient for vz(i-1,j-1,k)
                Lf(inp,invz - 4 - yn*4 - xn*yn*4, -kcont/zstp); % Coefficient for vz(i-1,j-1,k-1)
                
                % Right-hand-side part:0
                R(inp,1)                =   -a(i,j,k)*kcont;
            end
            
            % x-Stokes equation
            % Ghost vx unknowns (i=yn || k=zn) and boundary nodes (i=1, i=yn-1, j=1, j=xn, k=1, k=zn-1)
            if(i==yn || k==zn || i==1 || i==yn-1 || j==1 || j==xn || k==1 || k==zn-1)
                % Ghost vx unknowns (i=yn or k=zn: vx(i,j,k)=0
                if(i==yn || k==zn)
                    Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                    R(invx,1   )        =   0;          % Right-hand-side part
                end
                
                % Left and Right boundaries (j=1, j=xn)
                %             Left boundary
                if((j==1 ) && (i<yn) && (k<zn))
                    if (bc == DIRICHLET_BC)
                        %                 Dirichlet condition: vx(i,j,k) = 0
                        Lf(invx,invx, 1*kbond);  %Coeff. for vx(i,j,k)
                        R(invx,1)   = 0;
                    elseif (bc == NEUMANN_BC)
                        % Neumann condition: -vx(i,j+1,k)+vx(i,j,k)=0
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx+yn*4, -1*kbond);    % coeff. for vx(i,j+1,k)
                        R(invx,1)           =   0;          % Right-hand-side part
                        %
                        %                 elseif (bc == PSI_BC)
                        %                     %                 psi boudnary condition
                        %                     L(invx,invx) = 1*kbond;
                        %                     R(invx,1) = kbond*(psi(i,j)/xstp);
                    end
                end
                %             Right boundary
                if((j==xn) && (i<yn) && (k<zn))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition: vx(i,j,k) = 0
                        Lf(invx,invx, 1*kbond); %Coeff. for vx(i,j,k)
                        R(invx,1)    = 0;
                    elseif(bc==NEUMANN_BC)
                        % %                 Neumann condition: vx(i,j,k) - vx(i,j-1,k) = 0
                        Lf(invx,invx, 1*kbond);  %Coeff. for vx(i,j,k)
                        Lf(invx,invx-yn*4, -1*kbond); %Coeff. for vx(i,j-1,k)
                        R(invx,1) = 0;
                        %                 elseif(bc==PSI_BC)
                        %                     %             psi boundary condition
                        %                     L(invx,invx) = 1*kbond;
                        %                     R(invx,1) = kbond*(psi(i,j)/xstp);
                    end
                end
                %             Upper and lower boundaries
                % Upper boundary, inner points w.r.t x-axis, (i=1, k<zn, 1<j<xn)
                if(i==1 && k<zn && j>1 && j<xn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition: vx(i,j,k) = 0
                        % % equivalent to No slip vx=0: vx(i,j,k)-1/3*vx(i+1,j,k)=0.
                        %                     this extrapolation assumption might have affect on
                        %                     not matching the divergence theorem exactly!!
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx+4, -1/3*kbond); % Coefficient for vx(i+1,j,k)
                        R(invx,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i+1,j,k)-vx(i,j,k)=0
                        %                 Same for psi boundary condition!
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx+4, -1*kbond);    % Coefficient for vx(i+1,j,k)
                        R(invx,1)           =   0;          % Right-hand-side part
                    end
                end
                % Lower boundary, inner points w.r.t x-axis (i=yn-1, k<zn, 1<j<xn)
                if(i==yn-1 && k<zn && j>1 && j<xn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vx=0: vx(i,j,k)-1/3*vx(i-1,j,k)=0
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx-4, -1/3*kbond); % Coefficient for vx(i-1,j,k)
                        R(invx,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i,j)-vx(i-1,j)=0
                        %                 of psi condition
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx-4, -1*kbond); % Coefficient for vx(i-1,j,k)
                        R(invx,1)           =   0; % Right-hand-side part
                    end
                end
                
                
                %             Front and back boundaries
                % Front boundary, inner points w.r.t, (k=1, 1<j<xn, 1<i<yn-1)
                if(k==1 && j>1 && j<xn && i>1 && (i<(yn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition: vx(i,j,k) = 0
                        % % equivalent to No slip vx=0: vx(i,j,k)-1/3*vx(i,j,k+1)=0.
                        %                     this extrapolation assumption might have affect on
                        %                     not matching the divergence theorem exactly!!
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx+ xn*yn*4, -1/3*kbond); % Coefficient for vx(i,j,k+1)
                        R(invx,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i,j,k+1)-vx(i,j,k)=0
                        %                 Same for psi boundary condition!
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx + xn*yn*4, -1*kbond);    % Coefficient for vx(i,j,k+1)
                        R(invx,1)           =   0;          % Right-hand-side part
                    end
                end
                % Back boundary, inner points (k=zn-1, 1<j<xnum, 1<i<yn-1)
                if(k==zn-1 && j>1 && j<xn && i>1 && (i<(yn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vx=0: vx(i,j,k)-1/3*vx(i,j,k-1)=0
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx - xn*yn*4, -1/3*kbond); % Coefficient for vx(i,j,k-1)
                        R(invx,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i,j,k)-vx(i,j,k-1)=0
                        %                 of psi condition
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx - xn*yn*4, -1*kbond); % Coefficient for vx(i,j,k-1)
                        R(invx,1)           =   0; % Right-hand-side part
                    end
                end
                
                
                %Internal nodes: mu(d2vx/dx2+d2vx/dy2)-dP/dx=0
            else
                %dSxx/dx=2*muc(i+1,j+1,k+1)*(vx(i,j+1,k)-vx(i,j,k))/dx^2-2*muc(i+1,j,k+1)*(vx(i,j,k)-vx(i,j-1,k))/dx^2
                Lf(invx,invx + yn*4, 2*muc(i+1,j+1,k+1)/xstp^2);                     % Coefficient for vx(i,j+1,k)
                Lf(invx,invx - yn*4, 2*muc(i+1,j,k+1)/xstp^2);                           % Coefficient for vx(i,j-1,k)
                Lf(invx,invx, -2*muc(i+1,j+1,k+1)/xstp^2 - 2*muc(i+1,j,k+1)/xstp^2);   % Coefficient for vx(i,j,k)
                
                %dSxy/dy=muxy(i+1,j,k)*((vx(i+1,j,k)-vx(i,j,k))/dy^2+(vy(i+1,j,k)-vy(i+1,j-1,k))/dx/dy)-
                %         -muxy(i,j,k)*((vx(i,j,k)-vx(i-1,j,k))/dy^2+(vy(i,j,k)-vy(i,j-1,k))/dx/dy)-
                Lf(invx,invx+4, muxy(i+1,j,k)/ystp^2);                             % Coefficient for vx(i+1,j,k)
                Lf(invx,invx-4, muxy(i,j,k)/ystp^2);                               % Coefficient for vx(i-1,j,k)
                Lf(invx,invx, -muxy(i+1,j,k)/ystp^2-muxy(i,j,k)/ystp^2); % ADD coefficient for vx(i,j,k)
                Lf(invx,invy+4, muxy(i+1,j,k)/xstp/ystp);                          % Coefficient for vy(i+1,j,k)
                Lf(invx,invy+4-yn*4, -muxy(i+1,j,k)/xstp/ystp);                         % Coefficient for vy(i+1,j-1,k)
                Lf(invx,invy, -muxy(i,j,k)/xstp/ystp);                           % Coefficient for vy(i,j,k)
                Lf(invx,invy-yn*4, muxy(i,j,k)/xstp/ystp);                            % Coefficient for vy(i,j-1,k)
                
                %dSxz/dz=muxz(i,j,k+1)*((vx(i,j,k+1)-vx(i,j,k))/dz^2+(vz(i,j,k+1)-vz(i,j-1,k+1))/dx/dz)-
                %         -muxz(i,j,k)*((vx(i,j,k)-vx(i,j,k-1))/dz^2+(vz(i,j,k)-vz(i,j-1,k))/dx/dz)-
                Lf(invx,invx + xn*yn*4, muxz(i,j,k+1)/zstp^2);                             % Coefficient for vx(i,j,k+1)
                Lf(invx,invx - xn*yn*4, muxz(i,j,k)/zstp^2);                               % Coefficient for vx(i,j,k-1)
                Lf(invx,invx, -muxz(i,j,k+1)/zstp^2-muxz(i,j,k)/zstp^2); % ADD coefficient for vx(i,j,k)
                Lf(invx,invz+ xn*yn*4, muxz(i,j,k+1)/xstp/zstp);                          % Coefficient for vz(i,j,k+1)
                Lf(invx,invz - yn*4 + xn*yn*4, -muxz(i+1,j,k)/xstp/zstp);                         % Coefficient for vz(i,j-1,k+1)
                Lf(invx,invz, -muxz(i,j,k)/xstp/zstp);                           % Coefficient for vz(i,j,k)
                Lf(invx,invz - yn*4, muxz(i,j,k)/xstp/zstp);                            % Coefficient for vz(i,j-1,k)
                
                
                % -dP/dx=(P(i+1,j,k+1)-P(i+1,j+1,k+1))/dx
                Lf(invx,inp+4+xn*yn*4, kcont/xstp);                                     % Coefficient for P(i+1,j,k+1)
                Lf(invx,inp+4+yn*4+xn*yn*4, -kcont/xstp);                                    % Coefficient for P(i+1,j+1,k+1)
                % Right part:da/dx = (mu+lambda)*(a(i+1,j+1,k+1)-a(i+1,j,k+1))/dx
                R(invx,1)               =   ((muc(i+1,j+1,k+1)+muc(i+1,j,k+1))/2 ...
                    + (lambdac(i+1,j+1,k+1)+lambdac(i+1,j,k+1))/2) * (a(i+1,j+1,k+1)-a(i+1,j,k+1))/xstp;
            end
            
            % y-Stokes equation
            % Ghost vy unknowns (j=xn || k=zn) and boundary nodes (i=1, i=yn, j=1, j=xn-1, k=1, k=zn-1)
            if(j==xn || k==zn || i==1 || i==yn || j==1 || j==xn-1 || k==1 || k==zn-1)
                % Ghost vy unknowns (j=xn or k=zn: vy(i,j)=0
                if(j==xn || k==zn)
                    Lf(invy,invy, 1*kbond);                    % Coefficient for vy(i,j,k)
                    R(invy,1)           =   0;
                end
                % Upper boundary, i=1
                if((i==1) && (j<xn) && (k<zn))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition vy = 0
                        Lf(invy,invy, 1*kbond);
                        R(invy,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %            Neumann codition: -vy(i+1,j,k) + vy(i,j,k) = 0;
                        Lf(invy,invy, 1*kbond);                    % Coefficient for vy(i,j,k)
                        Lf(invy,invy+4, -1*kbond);    %Coefficient for vy(i+1,j,k)
                        R(invy,1)           =   0;
                        %                 elseif(bc==PSI_BC)
                        %                     %             psi boundary condtion:
                        %                     L(invy,invy) = 1*kbond;
                        %                     R(invy,1) = kbond*(psi(i,j)/ystp);
                    end
                end
                %             Lower boundary, i=yn
                if(i==yn && j<xn && k<zn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition vy = 0
                        Lf(invy,invy, 1*kbond);
                        R(invy,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %           Neumann codition: vy(i,j,k) - vy(i-1,j,k) = 0;
                        Lf(invy,invy, 1*kbond);   %Coeff. for vy(i,j,k)
                        Lf(invy, invy - 4, -1*kbond); %Coeff. for vy(i-1,j,k)
                        R(invy,1) = 0;               %Right side
                        %                 elseif(bc==PSI_BC)
                        %                     %                 psi boundary condtion:
                        %                     L(invy,invy) = 1*kbond;
                        %                     R(invy,1) = kbond*(psi(i,j)/ystp);
                    end
                end
                
                %             Left and right boundaries
                % Left boundary, inner points w.r.t. y-axis, (j=1, k<zn, 1<i<yn,)
                if(j==1 && k<zn && i>1 && i<yn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or eqv. to no slip:
                        %             % No slip vy=0: vy(i,j,k)-1/3*vy(i,j+1,k)=0
                        Lf(invy,invy, 1*kbond); % Coefficient for vy(i,j,k)
                        Lf(invy,invy+yn*4, -1/3*kbond); % Coefficient for vy(i,j+1,k)
                        R(invy,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                     elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann dvy/dx=0: vy(i,j,k)-vy(i,j+1,k)=0
                        %                 eqv. for psi boundary condtion
                        Lf(invy,invy, 1*kbond);                   % Coefficient for vy(i,j,k)
                        Lf(invy,invy+yn*4, -1*kbond);                   % Coefficient for vy(i,j+1,k)
                        R(invy,1)=0;
                    end
                end
                % Right boundary, inner points w.r.t y-axis, (j=xn-1, k<zn, 1<i<yn)
                if(j==xn-1 && k<zn && i>1 && i<yn)
                    if(bc==DIRICHLET_BC)
                        %                 Dirichlet or no slip: vy=0: vy(i,j,k)-1/3*vy(i,j-1,k)=0
                        Lf(invy,invy, 1*kbond); % Coefficient for vy(i,j,k)
                        Lf(invy,invy-yn*4, -1/3*kbond); % Coefficient for vy(i,j-1,k)
                        R(invy,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc==PSI_BC)
                        % Free slip dvy/dx=0: vy(i,j,k)-vy(i,j-1,k)=0
                        %                 eqv. to psi boundary condition
                        Lf(invy,invy, 1*kbond);                    % Coefficient for vy(i,j)
                        Lf(invy,invy-yn*4, -1*kbond);                   % Coefficient for vy(i,j-1)
                        R(invy,1)=0;
                    end
                end
                
                %             Front and back boundaries
                % Front boundary, inner points (k=1, 1<j<xn-1, 1<i<yn)
                if(k==1 && i>1 && i<yn && j>1 && (j<(xn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition: vy(i,j,k) = 0
                        % % equivalent to No slip vy=0: vy(i,j,k)-1/3*vy(i,j,k+1)=0.
                        %                     this extrapolation assumption might have affect on
                        %                     not matching the divergence theorem exactly!!
                        Lf(invy,invy, 1*kbond);    % Coefficient for vy(i,j,k)
                        Lf(invy,invy + xn*yn*4, -1/3*kbond); % Coefficient for vy(i,j,k+1)
                        R(invy,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvy/dz=0: vy(i,j,k)-vy(i,j,k+1)=0
                        %                 Same for psi boundary condition!
                        Lf(invy,invy, 1*kbond);    % Coefficient for vy(i,j,k)
                        Lf(invy,invy + xn*yn*4, -1*kbond);    % Coefficient for vy(i,j,k+1)
                        R(invy,1)           =   0;          % Right-hand-side part
                    end
                end
                % Back boundary, inner points (k=zn-1, 1<j<xn-1, 1<i<yn)
                if(k==zn-1 && i>1 && i<yn && j>1 && (j<(xn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vy=0: vy(i,j,k)-1/3*vy(i,j,k-1)=0
                        Lf(invy,invy, 1*kbond); % Coefficient for vy(i,j,k)
                        Lf(invy,invy - xn*yn*4, -1/3*kbond); % Coefficient for vy(i,j,k-1)
                        R(invy,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvy/dz=0: vy(i,j,k)-vy(i,j,k-1)=0
                        %                 of psi condition
                        Lf(invy,invy, 1*kbond); % Coefficient for vy(i,j,k)
                        Lf(invy,invy - xn*yn*4, -1*kbond); % Coefficient for vy(i,j,k-1)
                        R(invy,1)           =   0; % Right-hand-side part
                    end
                end
                
                
                %Internal nodes:
                %mu(d2vy/dx2+d2vy/dy2+d2vy/dz2)-dP/dy=(mu+lambda)da/dy
            else
                %dSyy/dy=2*muc(i+1,j+1,k+1)*(vy(i+1,j,k)-vy(i,j,k))/dy^2-2*muc(i,j+1,k+1)*(vy(i,j,k)-vy(i-1,j,k))/dy^2
                Lf(invy,invy+4, 2*muc(i+1,j+1,k+1)/ystp^2);                 % Coefficient for vy(i+1,j,k)
                Lf(invy,invy-4, 2*muc(i,j+1,k+1)/ystp^2);                   % Coefficient for vy(i-1,j,k)
                Lf(invy,invy, -2*muc(i+1,j+1,k+1)/ystp^2-2*muc(i,j+1,k+1)/ystp^2); % Coefficient for vy(i,j,k)
                
                %dSyx/dx=muyx(i,j+1,k)*((vy(i,j+1,k)-vy(i,j,k))/dx^2+(vx(i,j+1,k)-vx(i-1,j+1,k))/dx/dy)-
                %         -muyx(i,j,k)*((vy(i,j,k)-vy(i,j-1,k))/dx^2+(vx(i,j,k)-vx(i-1,j,k))/dx/dy)-
                Lf(invy,invy+yn*4, muyx(i,j+1,k)/xstp^2);                     % Coefficient for vy(i,j+1,k)
                Lf(invy,invy-yn*4, muyx(i,j,k)/xstp^2);                       % Coefficient for vy(i,j-1,k)
                Lf(invy,invy, -muyx(i,j+1,k)/xstp^2-muyx(i,j,k)/xstp^2); % ADD coefficient for vy(i,j,k)
                Lf(invy,invx+yn*4, muyx(i,j+1,k)/xstp/ystp);                  % Coefficient for vx(i,j+1,k)
                Lf(invy,invx-4+yn*4, -muyx(i,j+1,k)/xstp/ystp);                 % Coefficient for vx(i-1,j+1,k)
                Lf(invy,invx, -muyx(i,j,k)/xstp/ystp);                   % Coefficient for vx(i,j,k)
                Lf(invy,invx-4, muyx(i,j,k)/xstp/ystp);                    % Coefficient for vx(i-1,j,k)
                
                %dSyz/dz=muyz(i,j,k+1)*((vy(i,j,k+1)-vy(i,j,k))/dz^2+(vz(i,j,k+1)-vz(i-1,j,k+1))/dz/dy)-
                %         -muyz(i,j,k)*((vy(i,j,k)-vy(i,j,k-1))/dz^2+(vz(i,j,k)-vz(i-1,j,k))/dz/dy)-
                Lf(invy,invy+xn*yn*4, muyz(i,j,k+1)/zstp^2);                     % Coefficient for vy(i,j,k+1)
                Lf(invy,invy-xn*yn*4, muyz(i,j,k)/zstp^2);                       % Coefficient for vy(i,j,k-1)
                Lf(invy,invy, -muyz(i,j,k+1)/zstp^2-muyz(i,j,k)/zstp^2); % ADD coefficient for vy(i,j,k)
                Lf(invy,invz+xn*yn*4, muyz(i,j,k+1)/zstp/ystp);                  % Coefficient for vz(i,j,k+1)
                Lf(invy,invz-4+xn*yn*4, -muyz(i,j,k+1)/zstp/ystp);                 % Coefficient for vz(i-1,j,k+1)
                Lf(invy,invz, -muyz(i,j,k)/zstp/ystp);                   % Coefficient for vz(i,j,k)
                Lf(invy,invz-4, muyz(i,j,k)/zstp/ystp);                    % Coefficient for vz(i-1,j,k)
                
                % -dP/dy=(P(i,j+1,k+1)-P(i+1,j+1,k+1))/dy
                Lf(invy,inp+yn*4+xn*yn*4, kcont/ystp);                             % Coefficient for P(i,j+1,k+1)
                Lf(invy,inp+4+yn*4+xn*yn*4, -kcont/ystp);                            % Coefficient for P(i+1,j+1,k+1)
                % Right part: da/dy =
                % (mu+lambda)*(a(i+1,j+1,k+1)-a(i,j+1,k+1))/dy
                R(invy,1)               = ((muc(i+1,j+1,k+1)+muc(i,j+1,k+1))/2 ...
                    + (lambdac(i+1,j+1,k+1)+lambdac(i,j+1,k+1))/2)*(a(i+1,j+1,k+1)-a(i,j+1,k+1))/ystp;
            end
            
            
            % z-Stokes equation
            % Ghost vz unknowns (i=yn || j=xn) and boundary nodes (i=1, i=yn-1, j=1, j=xn-1, k=1, k=zn)
            if(i==yn || j==xn ||  k==1 || k==zn || j==1 || j==xn-1 || i==1 || i==yn-1)
                % Ghost vz unknowns (j=xn or i=yn: vz(i,j,k)=0
                if( i==yn || j==xn)
                    Lf(invz,invz, 1*kbond);                    % Coefficient for vz(i,j,k)
                    R(invz,1)           =   0;
                end
                % Front boundary
                if((k==1) && (i<yn) && (j<xn))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition vz = 0
                        Lf(invz,invz, 1*kbond);
                        R(invz,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %            Neumann codition: -vz(i,j,k+1) + vy(i,j,k) = 0;
                        Lf(invz,invz, 1*kbond);                    % Coefficient for vz(i,j,k)
                        Lf(invz,invz+xn*yn*4, -1*kbond);    %Coefficient for vz(i,j,k+1)
                        R(invz,1)           =   0;
                        %                 elseif(bc==PSI_BC)
                        %                     %             psi boundary condtion:
                        %                     L(invy,invy) = 1*kbond;
                        %                     R(invy,1) = kbond*(psi(i,j)/ystp);
                    end
                end
                %  Back boundary, k=zn
                if(k==zn && i<yn && j<xn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition vz = 0
                        Lf(invz,invz, 1*kbond);
                        R(invz,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %           Neumann codition: vz(i,j,k) - vz(i,j,k-1) = 0;
                        Lf(invz,invz, 1*kbond);   %Coeff. for vz(i,j,k)
                        Lf(invz, invz - xn*yn*4, -1*kbond); %Coeff. for vz(i,j,k-1)
                        R(invz,1) = 0;               %Right side
                    end
                end
                
                %             Left and right boundaries
                % Left boundary, inner points w.r.t. z-axis, (j=1, i<yn, 1<k<zn,)
                if(j==1 && i<yn && k>1 && k<zn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or eqv. to no slip:
                        %             % No slip vz=0: vz(i,j,k)-1/3*vz(i,j+1,k)=0
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz+yn*4, -1/3*kbond); % Coefficient for vz(i,j+1,k)
                        R(invz,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                     elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann dvz/dx=0: vz(i,j,k)-vz(i,j+1,k)=0
                        %                 eqv. for psi boundary condtion
                        Lf(invz,invz, 1*kbond);                   % Coefficient for vz(i,j,k)
                        Lf(invz,invz+yn*4, -1*kbond);                   % Coefficient for vz(i,j+1,k)
                        R(invz,1)=0;
                    end
                end
                % Right boundary, inner points w.r.t z-axis, (j=xn-1, i<yn, 1<k<zn)
                if(j==xn-1 && i<yn && k>1 && k<zn)
                    if(bc==DIRICHLET_BC)
                        %                 Dirichlet or no slip: vz=0: vz(i,j,k)-1/3*vz(i,j-1,k)=0
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz-yn*4, -1/3*kbond); % Coefficient for vz(i,j-1,k)
                        R(invz,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc==PSI_BC)
                        % Free slip dvz/dx=0: vz(i,j,k)-vz(i,j-1,k)=0
                        %                 eqv. to psi boundary condition
                        Lf(invz,invz, 1*kbond);                    % Coefficient for vz(i,j,k)
                        Lf(invz,invz-yn*4, -1*kbond);                   % Coefficient for vz(i,j-1,k)
                        R(invz,1)=0;
                    end
                end
                
                %   Upper and Lower Boundaries:
                % Upper boundary, inner points (i=1, 1<j<xn-1, 1<k<zn)
                if(i==1 && k>1 && k<zn && j>1 && (j<(xn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition: vz(i,j,k) = 0
                        % % equivalent to No slip vz=0: vz(i,j,k)-1/3*vz(i+1,j,k)=0.
                        %                     this extrapolation assumption might have affect on
                        %                     not matching the divergence theorem exactly!!
                        Lf(invz,invz, 1*kbond);    % Coefficient for vz(i,j,k)
                        Lf(invz,invz + 4, -1/3*kbond); % Coefficient for vy(i+1,j,k)
                        R(invz,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvz/dy=0: vz(i+1,j,k)-vz(i,j,k)=0
                        %                 Same for psi boundary condition!
                        Lf(invz,invz, 1*kbond);    % Coefficient for vz(i,j,k)
                        Lf(invz,invz + 4, -1*kbond);    % Coefficient for vz(i+1,j,k)
                        R(invz,1)           =   0;          % Right-hand-side part
                    end
                end
                % Lower boundary, inner points (i=yn-1, 1<j<xn-1, 1<k<zn)
                if((i==(yn-1)) && k>1 && k<zn && j>1 && (j<(xn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vz=0: vz(i,j,k)-1/3*vz(i-1,j,k)=0
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz - 4, -1/3*kbond); % Coefficient for vz(i-1,j,k)
                        R(invz,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvz/dy=0: vz(i,j,k)-vz(i-1,j,k)=0
                        %                 of psi condition
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz - 4, -1*kbond); % Coefficient for vz(i-1,j,k)
                        R(invz,1)           =   0; % Right-hand-side part
                    end
                end
                
                %Internal nodes:
                %mu(d2vz/dx2+d2vz/dy2+d2vz/dz2)-dP/dz=(mu+lambda)da/dz
            else
                %dSzz/dz=2*muc(i+1,j+1,k+1)*(vz(i,j,k+1)-vz(i,j,k))/dz^2-2*muc(i+1,j+1,k)*(vz(i,j,k)-vz(i,j,k-1))/dz^2
                Lf(invz,invz+xn*yn*4, 2*muc(i+1,j+1,k+1)/zstp^2);                 % Coefficient for vz(i,j,k+1)
                Lf(invz,invz-xn*yn*4, 2*muc(i+1,j+1,k)/zstp^2);                   % Coefficient for vz(i,j,k-1)
                Lf(invz,invz, -2*muc(i+1,j+1,k+1)/zstp^2-2*muc(i+1,j+1,k)/zstp^2); % Coefficient for vz(i,j,k)
                
                %dSzx/dx=muzx(i,j+1,k)*((vz(i,j+1,k)-vz(i,j,k))/dx^2+(vx(i,j+1,k)-vx(i,j+1,k-1))/dx/dz)-
                %         -muzx(i,j,k)*((vz(i,j,k)-vz(i,j-1,k))/dx^2+(vx(i,j,k)-vx(i,j,k-1))/dx/dz)-
                Lf(invz,invz+yn*4, muzx(i,j+1,k)/xstp^2);                     % Coefficient for vz(i,j+1,k)
                Lf(invz,invz-yn*4, muzx(i,j,k)/xstp^2);                       % Coefficient for vz(i,j-1,k)
                Lf(invz,invz, -muzx(i,j+1,k)/xstp^2-muzx(i,j,k)/xstp^2); % ADD coefficient for vz(i,j,k)
                Lf(invz,invx+yn*4, muzx(i,j+1,k)/xstp/zstp);                  % Coefficient for vx(i,j+1,k)
                Lf(invz,invx+yn*4-xn*yn*4, -muzx(i,j+1,k)/xstp/zstp);                 % Coefficient for vx(i,j+1,k-1)
                Lf(invz,invx, -muzx(i,j,k)/xstp/zstp);                   % Coefficient for vx(i,j,k)
                Lf(invz,invx-xn*yn*4, muzx(i,j,k)/xstp/zstp);                    % Coefficient for vx(i,j,k-1)
                
                %dSzy/dz=muyz(i+1,j,k)*((vz(i+1,j,k)-vz(i,j,k))/dy^2+(vy(i+1,j,k)-vy(i+1,j,k-1))/dy/dz)-
                %         -muyz(i,j,k)*((vz(i,j,k)-vz(i-1,j,k))/dy^2+(vy(i,j,k)-vy(i,j,k-1))/dy/dz)-
                Lf(invz,invz+4, muyz(i+1,j,k)/ystp^2);                     % Coefficient for vz(i+1,j,k)
                Lf(invz,invz-4, muyz(i,j,k)/ystp^2);                       % Coefficient for vz(i-1,j,k)
                Lf(invz,invz, -muyz(i+1,j,k)/ystp^2-muyz(i,j,k)/ystp^2); % ADD coefficient for vz(i,j,k)
                Lf(invz,invy+4, muyz(i+1,j,k)/zstp/ystp);                  % Coefficient for vy(i+1,j,k)
                Lf(invz,invy+4-xn*yn*4, -muyz(i+1,j,k)/zstp/ystp);                 % Coefficient for vy(i+1,j,k-1)
                Lf(invz,invy, -muyz(i,j,k)/zstp/ystp);                   % Coefficient for vy(i,j,k)
                Lf(invz,invy-xn*yn*4, muyz(i,j,k)/zstp/ystp);                    % Coefficient for vy(i,j,k-1)
                
                % -dP/dz=(P(i+1,j+1,k)-P(i+1,j+1,k+1))/dz
                Lf(invz,inp+4+yn*4, kcont/zstp);                             % Coefficient for P(i+1,j+1,k)
                Lf(invz,inp+4+yn*4+xn*yn*4, -kcont/zstp);                            % Coefficient for P(i+1,j+1,k+1)
                % Right part: da/dz =
                % (mu+lambda)*(a(i+1,j+1,k+1)-a(i+1,j+1,k))/dy
                R(invz,1)          = ((muc(i+1,j+1,k+1)+muc(i+1,j+1,k))/2 ...
                    + (lambdac(i+1,j+1,k+1)+lambdac(i+1,j+1,k))/2)*(a(i+1,j+1,k+1)-a(i+1,j+1,k))/zstp;
            end
        end
    end
end



% Constraints when neumann boundary for displacement is used for all the boundaries
if (bc==NEUMANN_BC)
    % For rigid body translation
    % % Let's set the average of the each of the displacement components to be
    % % zero!
    for i = 2:yn-1
        for j = 2:xn-1
            for k = 2:zn-1
                invx = ((k-1)*xn*yn + (j-1)*yn + i)*4 - 2;
                invy = invx + 1;
                invz = invx + 2;
                Lf(1+NS,invx, 1*kbond); % Coefficient for vx(i,j,k)
                Lf(2+NS,invy, 1*kbond); % Coefficient for vy(i,j,k)
                Lf(3+NS,invz, 1*kbond); % Coefficient for vz(i,j,k)
            end
        end
    end
    R(1+xn*yn*4,1) = 0;
    R(2+xn*yn*4,1) = 0;
    R(3+xn*yn*4,1) = 0;
    
    
    % For rigid body rotation
    % At all interior points, mean of rigid body rotation tensor components set
    % to zero:
    % i.e. Omegaxy = 0 => dvx/dy - dvy/dx = 0;
    %     Omegaxz = 0 => dvx/dz - dvz/dx = 0;
    %       Omegayz = 0 => dvy/dz - dvz/dy = 0;
    % Now considering Omega lying in the grid in the same way as S (sigma),
    % For Oxy: (vx(i,j,k)-vx(i-1,j,k))/ystp - (vy(i,j,k)-vy(i,j-1,k))/xstp = 0;
    % For Oxz: (vx(i,j,k)-vx(i,j,k-1))/zstp - (vz(i,j,k)-vz(i,j-1,k))/xstp = 0;
    % For Oyz: (vy(i,j,k)-vy(i,j,k-1))/zstp - (vz(i,j,k)-vz(i-1,j,k))/ystp = 0;
    for i = 3:yn-1
        for j = 3:xn-1
            for k = 3:zn-1
                invx = ((k-1)*xn*yn + (j-1)*yn+i)*4 - 2;
                invy = invx + 1;
                invz = invx + 2;
                %             For Oxy
                Lf(4+NS,invx, kbond/ystp); % Coefficient for vx(i,j,k)
                Lf(4+NS,invx-4, -kbond/ystp); % Coefficient for vx(i-1,j,k)
                Lf(4+NS,invy, -kbond/xstp); % Coefficient for vy(i,j,k)
                Lf(4+NS,invy-yn*4, kbond/xstp); % Coefficient for vy(i,j-1,k)
                
                %             For Oxz
                Lf(5+NS,invx, kbond/zstp); % Coefficient for vx(i,j,k)
                Lf(5+NS,invx-xn*yn*4, -kbond/zstp); % Coefficient for vx(i,j,k-1)
                Lf(5+NS,invz, -kbond/xstp); % Coefficient for vz(i,j,k)
                Lf(5+NS,invz-yn*4, kbond/xstp); % Coefficient for vz(i,j-1,k)
                
                %             For Oxy
                Lf(6+NS,invy, kbond/zstp); % Coefficient for vy(i,j,k)
                Lf(6+NS,invy-xn*yn*4, -kbond/zstp); % Coefficient for vy(i,j,k-1)
                Lf(6+NS,invz, -kbond/ystp); % Coefficient for vz(i,j,k)
                Lf(6+NS,invz-4, kbond/ystp); % Coefficient for vz(i-1,j,k)
                
            end
        end
    end
    R(4+NS,1) = 0;
    R(5+NS,1) = 0;
    R(6+NS,1) = 0;
    
    ecnstrnt_indx = 7;
    
    if (pChoice == ENFORCE_CC) %Since CC is enforced, extra constraint needed here to fix pressure value to make it full rank!
        inp     =   ((p0cellz-1)*xn*yn + (p0cellx-1)*yn + p0celly)*4 - 3;
        Lf(ecnstrnt_indx + NS,inp, 1*kbond);  %Coeff. for p0cell
        R(ecnstrnt_indx +NS,1) = p0cell;
        ecnstrnt_indx = ecnstrnt_indx + 1;
    end
    
    
end

if (bc ~= NEUMANN_BC && pChoice == ENFORCE_CC)
    inp     =   ((p0cellz-1)*xn*yn + (p0cellx-1)*yn + p0celly)*4 - 3;
    Lf(1+NS,inp, 1*kbond);  %Coeff. for p0cell
    R(1+NS,1) = p0cell;
    ecnstrnt_indx = ecnstrnt_indx + 1;
end

% Create sparse matrix from the indices and values
L = sparse(i_indx(1:vals_cntr),j_indx(1:vals_cntr),L_vals(1:vals_cntr), ecnstrnt+NS, NS);

toc;
display('matrix built');
tic;
%Obtaining vector of solutions S()
S=L\R;

% tol = 1e-12; maxit = 40;
% % setup.type = 'ilutp';
% % setup.udiag = 1;
% % [LL, UU] = ilu(L, setup);
% % [LL,UU] = ilu(L,struct('type','ilutp','droptol',1e-5,'milu','off','udiag',1));
% [LL,UU] = ilu(L,struct('type','ilutp','droptol',1e-5));
% [S,fl1,rr1,it1,rv1] = bicgstab(L,R,tol,maxit,LL,UU);
% % [S, fl1, rr1, it1, rv1] = bicgstab(L,R);
toc;
display('system solver timed');


% Reload solutions to 3D p(), vx(), vy(), vz() arrays
% Dimensions of arrays are reduced compared to the basic grid
p=zeros(yn,xn,zn);
vy=zeros(yn,xn,zn);
vx=zeros(yn,xn,zn);
vz=zeros(yn,xn,zn);
% Process all Grid points
for i=1:1:yn
    for j=1:1:xn
        for k = 1:1:zn
            % Global index for P, vx, vy and vz in S()
            inp=((k-1)*xn*yn + (j-1)*yn+i)*4 - 3; % P
            invx=inp+1;
            invy=inp+2;
            invz=inp+3;
            % P
            p(i,j,k)=S(inp)*kcont;
            % vx
            vx(i,j,k)=S(invx);
            % vy
            vy(i,j,k)=S(invy);
            %         vz
            vz(i,j,k)=S(invz);
        end
    end
end

[ax ay az] = meshgrid(1:xn,1:yn,1:zn);

% savevtk(p,'pressure3D.vtk');

%
% Compute vx,vy for internal nodes
vx1=zeros(yn,xn,zn);
vy1=zeros(yn,xn,zn);
vz1=zeros(yn,xn,zn);
% Process internal Grid points
for i=2:1:yn-1
    for j=2:1:xn-1
        for k = 2:1:zn-1
            % vx
            vx1(i,j,k)=(vx(i-1,j,k)+vx(i,j,k))/2;
            % vy
            vy1(i,j,k)=(vy(i,j-1,k)+vy(i,j,k))/2;
            % vz
            vz1(i,j,k)=(vz(i,j,k-1)+vz(i,j,k))/2;
        end
    end
end

vec3DToVtk(vx,vy,vz,ax,ay,az,'varmu3D.vtk');
vec3DToVtk(vx1,vy1,vz1,ax,ay,az,'varmu3D1.vtk');
save('results.mat','vx1','vy1','vz1','p', 'L','R');

end
