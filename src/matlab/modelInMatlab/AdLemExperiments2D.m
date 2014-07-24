% Solution of the Brain deformation with irregular domain (or brain-CSF
% considered as an interface and extending the boundary to the skull.
% Compressibility constraint enforced only on the brain region.
% Regular staggered grid using a pressure-velocity formulation
% for a medium with varying lame's parameters.
%
% Always expetcs a seg_mask with at least 1 width of CSF in the periphery.
% extra_CSF is the width of CSF around the provided mask.
% VENTRICLES_UNIFORM will adapt the atrophy data in div_slice by padding
% uniform atrophy in the periphery
function [v2 p] = AdLemExperiments2D(seg_mask,div_slice, extra_CSF, ...
    muBrain,muRatio,lambdaBrain,lambdaRatio, VENTRICLES_UNIFORM)


% Since mu and lambda Brain are greater than CSF, the ratio muB/muC given
% as input in the form muRatio should be greater or equal to 1.
if(lambdaRatio < 1 || muRatio < 1)
    error('invalid muRatio or invalid lamdaRatio, they should be greater than one');
end

[brain_mask min_max_row_col] = getRectConvexHull(seg_mask);
[ynum xnum] = size(brain_mask);

a = div_slice(min_max_row_col(1) : min_max_row_col(2), min_max_row_col(3) : min_max_row_col(4));
%pad with extra_CSF width of zeros all around.
brain_mask = padarray(brain_mask,[extra_CSF extra_CSF]);

if(VENTRICLES_UNIFORM == false)
    % distributed source outside convex hull(padded areas), excluding 4
    % corners where CC is not enforced! Total pixels of padded areas of
    % width w is (xn+2w)(yn+2w)-xn.yn
    a = padarray(a,[extra_CSF extra_CSF],-sum(a(:))/(2*extra_CSF*(xnum + ynum + 2*extra_CSF)));
else
    a = padarray(a,[extra_CSF extra_CSF]);
    a(~brain_mask) = 0;
    %         a(brain_mask==0) = 0;
    brain_mask_n = ~brain_mask;
    a(~brain_mask) = -sum(a(:))/sum(brain_mask_n(:));
end
%Grid nodes have one dimension greater than the cell centers!
ynum = ynum+extra_CSF*2+1;
xnum = xnum+extra_CSF*2+1;


% add a dummy atrophy corresponding to dummy pressure values!
a = [zeros(1,xnum); zeros(ynum-1,1) a];

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
xsize   =   (xnum-1)*xres;        % Horizontal
ysize   =   (ynum-1)*yres;        % Vertical

% Grid step
xstp    =   xsize/(xnum-1); % Horizontal
ystp    =   ysize/(ynum-1); % Vertical

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
muc = zeros(ynum-1,xnum-1);
lambdac = zeros(ynum-1,xnum-1);

muc(brain_mask==true) = muBrain;
muc(brain_mask==false) = muCSF;
lambdac(brain_mask==true) = lambdaBrain;
lambdac(brain_mask==false) = lambdaCSF;

% Interpolate mu and lambda on basic nodes from values at centers of the cells
mug = (muc(:,1:xnum-2) + muc(:,2:xnum-1))/2; %horizontal interpolation
mug = (mug(1:ynum-2,:) + mug(2:ynum-1,:))/2; %vertical interpolation

% Pad the values on the boundary edges with that of adjacent ones.
% horizontal padding
mug = [mug(1,:); mug; mug(end,:)];
% vertical padding
mug = [mug(:,1) mug mug(:,end)];
% Now Add dummy zeros on 1st row and 1st col of centered cells:
muc = [zeros(1,xnum); zeros(ynum-1,1) muc];



lambdag = (lambdac(:,1:xnum-2) + lambdac(:,2:xnum-1))/2; %horizontal interpolation
lambdag = (lambdag(1:ynum-2,:) + lambdag(2:ynum-1,:))/2; %vertical interpolation

% Pad the values on the boundary edges with that of adjacent ones.
% horizontal padding
lambdag = [lambdag(1,:); lambdag; lambdag(end,:)];
% vertical padding
lambdag = [lambdag(:,1) lambdag lambdag(:,end)];

% Now Add dummy zeros on 1st row and 1st col of centered cells:
lambdac = [zeros(1,xnum); zeros(ynum-1,1) lambdac];


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
ecnstrnt_indx = 1;


% Top left:
p0cellx = 3; %indx j
p0celly = 2; %indx i

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
    psi = psi(2:ynum+1,2:xnum+1);
end


% Matrix of coefficients initialization
% L       =   sparse(ecnstrnt+xnum*ynum*3,xnum*ynum*3);
% Tentative number of non-zeros in L matrix
tot_non_zeros = 4*(xnum+ynum) + 4*(xnum-1)*(ynum-1) + 12*((ynum-3)*(xnum-2) + ...
    (xnum-3)*(ynum-2)) + 5 + 2*ynum + 2*xnum; %last ko 2*yn+2*xn ettikai extra
i_indx = zeros(tot_non_zeros,1);
j_indx = zeros(tot_non_zeros,1);
L_vals = zeros(tot_non_zeros,1);


% Vector of right part initialization
R       =   zeros(ecnstrnt+ xnum*ynum*3,1);


% Solving x-Stokes, y-Stokes and continuity equations
% x-Stokes: mu(d2vx/dx2+d2vx/dy2)-dP/dx=(mu+lambda)da/dx
% y-Stokes: mu(d2vy/dx2+d2vy/dy2)-dP/dy=gy*RHO + (mu+lambda)da/dy
% continuity: dvx/dx+dvy/dy=-a
% Compose matrix of coefficients L()
% and vector (column) of right parts R()
% Boundary conditions: free slip
% Process all Grid points

pnodes_cntr = 0;
vals_cntr = 0;
    function Lf(ii,jj,xx)
        vals_cntr = vals_cntr + 1;
        i_indx(vals_cntr) = ii;
        j_indx(vals_cntr) = jj;
        L_vals(vals_cntr) = xx;
    end

tic;
for i=1:1:ynum
    for j=1:1:xnum
        
        % Global index for P, vx, vy in the current node
        inp     =   ((j-1)*ynum+i)*3-2; % P
        invx    =   inp+1;
        invy    =   inp+2;
        
        %         It's the same for atrophy too, so we can use this index to access
        %         atrophy values too!
        
        % Continuity equation
%         True only on the brain region. We set ghost region and CSF region
%         to have p=0 and not enforcing the continuity equation. This
%         automatically takes care of four corners now since we have at
%         least two pixel width of CSF all around.

        % Ghost pressure unknowns (i=1, j=1) and boundary nodes (4 corners + one cell)
        if ( (i==1) || (j==1))
            % Ghost pressure unknowns (i=1, j=1): P(i,j)=0
                Lf(inp,inp,1*kbond);
                R(inp,1  )          =   0;          % Right-hand-side part
            %Internal nodes: dvx/dx+dvy/dy=-a
        else
%             CSF region set p=0, don't care about the continuity equation:
            if (brain_mask(i-1,j-1)==0)
                Lf(inp,inp,1*kbond);
                R(inp,1  )          =   0;          % Right-hand-side part
            else
                
            %dvx/dx=(vx(i-1,j)-vx(i-1,j-1))/dx
            Lf(inp,invx-3,kcont/xstp); % Coefficient for vx(i-1,j)
            Lf(inp,invx-3-ynum*3,-kcont/xstp); % Coefficient for vx(i-1,j-1)
            
            %dvy/dy=(vy(i,j-1)-vy(i-1,j-1))/dy
            Lf(inp,invy-ynum*3,kcont/ystp); % Coefficient for vy(i,j-1)
            Lf(inp,invy-3-ynum*3,-kcont/ystp); % Coefficient for vy(i-1,j-1)
            % Right-hand-side part:0
            R(inp,1)                =   -a(i,j)*kcont;
            pnodes_cntr = pnodes_cntr+1;
            end
        end
        
        % x-Stokes equation
        % Ghost vx unknowns (i=ynum) and boundary nodes (i=1, i=ynum-1, j=1, j=xnum)
        if(i==1 || i==ynum-1 || i==ynum || j==1 || j==xnum)
            % Ghost vx unknowns (i=ynum: vx(i,j)=0
            if(i==ynum)
                Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                R(invx,1   )        =   0;          % Right-hand-side part
            end
            % Left and Right boundaries (j=1, j=xnum)
            %             Left boundary
            if((j==1 ) && i<ynum)
                if (bc == DIRICHLET_BC)
                    %                 Dirichlet condition: vx(i,j) = 0
                    Lf(invx,invx,1*kbond);  %Coeff. for vx(i,j)
                    R(invx,1)   = 0;
                elseif (bc == NEUMANN_BC)
                    % Neumann condition: -vx(i,j+1)+vx(i,j)=0
                    Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                    Lf(invx,invx+ynum*3,-1*kbond);    % coeff. for vx(i,j+1)
                    R(invx,1)           =   0;          % Right-hand-side part
                    %
                elseif (bc == PSI_BC)
                    %                 psi boudnary condition
                    Lf(invx,invx,1*kbond);
                    R(invx,1) = kbond*(psi(i,j)/xstp);
                end
            end
            %             Right boundary
            if((j==xnum) && i<ynum)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition: vx(i,j) = 0
                    Lf(invx,invx,1*kbond); %Coeff. for vx(i,j)
                    R(invx,1)    = 0;
                elseif(bc==NEUMANN_BC)
                    % %                 Neumann condition: vx(i,j) - vx(i,j-1) = 0
                    Lf(invx,invx,1*kbond);  %Coeff. for vx(i,j)
                    Lf(invx,invx-ynum*3,-1*kbond); %Coeff. for vx(i,j-1)
                    R(invx,1) = 0;
                elseif(bc==PSI_BC)
                    %             psi boundary condition
                    Lf(invx,invx,1*kbond);
                    R(invx,1) = kbond*(psi(i,j)/xstp);
                end
            end
            % Upper boundary, inner points (i=1, 1<j<xnum)
            if(i==1 && j>1 && j<xnum)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition: vx(i,j) = 0
                    % % equivalent to No slip vx=0: vx(i,j)-1/3*vx(i+1,j)=0.
                    Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                    Lf(invx,invx+3,-1/3*kbond); % Coefficient for vx(i+1,j)
                    R(invx,1)          =   0;          % Right-hand-side part
                elseif(bc==NEUMANN_BC || bc == PSI_BC)
                    % Neumann condition dvx/dy=0: vx(i+1,j)-vx(i,j)=0
                    %                 Same for psi boundary condition!
                    Lf(invx,invx,1*kbond);    % Coefficient for vx(i,j)
                    Lf(invx,invx+3,-1*kbond);    % Coefficient for vx(i+1,j)
                    R(invx,1)           =   0;          % Right-hand-side part
                end
            end
            % Lower boundary, inner points (i=ynum-1, 1<j<xnum)
            if(i==ynum-1 && j>1 && j<xnum)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet or No slip:
                    % % No slip vx=0: vx(i,j)-1/3*vx(i-1,j)=0
                    Lf(invx,invx,1*kbond); % Coefficient for vx(i,j)
                    Lf(invx,invx-3,-1/3*kbond); % Coefficient for vx(i-1,j)
                    R(invx,1)         =   0; % Right part
                elseif(bc==NEUMANN_BC || bc == PSI_BC)
                    % Neumann condition dvx/dy=0: vx(i,j)-vx(i-1,j)=0
                    %                 of psi condition
                    Lf(invx,invx,1*kbond); % Coefficient for vx(i,j)
                    Lf(invx,invx-3,-1*kbond); % Coefficient for vx(i-1,j)
                    R(invx,1)           =   0; % Right-hand-side part
                end
            end
            %Internal nodes: mu(d2vx/dx2+d2vx/dy2)-dP/dx=0
        else
            %dSxx/dx=2*muc(i+1,j+1)*(vx(i,j+1)-vx(i,j))/dx^2-2*muc(i+1,j)*(vx(i,j)-vx(i,j-1))/dx^2
            Lf(invx,invx+ynum*3,2*muc(i+1,j+1)/xstp^2);                         % Coefficient for vx(i,j+1)
            Lf(invx,invx-ynum*3,2*muc(i+1,j)/xstp^2);                           % Coefficient for vx(i,j-1)
            Lf(invx,invx,-2*muc(i+1,j+1)/xstp^2-2*muc(i+1,j)/xstp^2);   % Coefficient for vx(i,j)
            
            %dSxy/dy=mug(i+1,j)*((vx(i+1,j)-vx(i,j))/dy^2+(vy(i+1,j)-vy(i+1,j-1))/dx/dy)-
            %         -mug(i,j)*((vx(i,j)-vx(i-1,j))/dy^2+(vy(i,j)-vy(i,j-1))/dx/dy)-
            Lf(invx, invx+3, mug(i+1,j)/ystp^2);                             % Coefficient for vx(i+1,j)
            Lf(invx, invx-3, mug(i,j)/ystp^2);                               % Coefficient for vx(i-1,j)
            Lf(invx, invx, -mug(i+1,j)/ystp^2-mug(i,j)/ystp^2); % ADD coefficient for vx(i,j)
            Lf(invx, invy+3, mug(i+1,j)/xstp/ystp);                          % Coefficient for vy(i+1,j)
            Lf(invx, invy+3-ynum*3, -mug(i+1,j)/xstp/ystp);                         % Coefficient for vy(i+1,j-1)
            Lf(invx, invy, -mug(i,j)/xstp/ystp);                           % Coefficient for vy(i,j)
            Lf(invx, invy-ynum*3, mug(i,j)/xstp/ystp);                            % Coefficient for vy(i,j-1)
            % -dP/dx=(P(i+1,j)-P(i+1,j+1))/dx
            Lf(invx, inp+3, kcont/xstp);                                     % Coefficient for P(i+1,j)
            Lf(invx, inp+3+ynum*3, -kcont/xstp);                                    % Coefficient for P(i+1,j+1)
            % Right part:da/dx = (a(i+1,j+1)-a(i+1,j))/dx
            R(invx,1)               =   ((muc(i+1,j+1)+muc(i+1,j))/2 + (lambdac(i+1,j+1)+lambdac(i+1,j))/2) * (a(i+1,j+1)-a(i+1,j))/xstp;
            
        end
        
        % y-Stokes equation
        % Ghost vy unknowns (j=xnum) and boundary nodes (i=1, i=ynum, j=1, j=xnum-1)
        if(i==1 || i==ynum || j==1 || j==xnum-1 || j==xnum)
            % Ghost vy unknowns (j=xnum: vy(i,j)=0
            if(j==xnum)
                Lf(invy,invy,1*kbond);                    % Coefficient for vy(i,j)
                R(invy,1)           =   0;
            end
            % Upper boundary, i=1
            if((i==1) && (j<xnum))
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition vy = 0
                    Lf(invy,invy,1*kbond);
                    R(invy,1) = 0;
                elseif(bc==NEUMANN_BC)
                    % %            Neumann codition: -vy(i+1,j) + vy(i,j) = 0;
                    Lf(invy,invy,1*kbond);                    % Coefficient for vy(i,j)
                    Lf(invy,invy+3,-1*kbond);    %Coefficient for vy(i+1,j)
                    R(invy,1)           =   0;
                elseif(bc==PSI_BC)
                    %             psi boundary condtion:
                    Lf(invy,invy,1*kbond);
                    R(invy,1) = kbond*(psi(i,j)/ystp);
                end
            end
            %             Lower boundary, i=ynum
            if(i==ynum && j < xnum)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet condition vy = 0
                    Lf(invy,invy,1*kbond);
                    R(invy,1) = 0;
                elseif(bc==NEUMANN_BC)
                    % %           Neumann codition: vy(i,j) - vy(i-1,j) = 0;
                    Lf(invy,invy,1*kbond);   %Coeff. for vy(i,j)
                    Lf(invy, invy - 3,-1*kbond); %Coeff. for vy(i-1,j)
                    R(invy,1) = 0;               %Right side
                elseif(bc==PSI_BC)
                    %                 psi boundary condtion:
                    Lf(invy,invy,1*kbond);
                    R(invy,1) = kbond*(psi(i,j)/ystp);
                end
            end
            
            % Left boundary, inner points (j=1, 1<i<ynum)
            if(j==1 && i>1 && i<ynum)
                if(bc==DIRICHLET_BC)
                    %                     Dirichlet or eqv. to no slip:
                    %             % No slip vy=0: vy(i,j)-1/3*vy(i,j+1)=0
                    Lf(invy,invy,1*kbond); % Coefficient for vy(i,j)
                    Lf(invy,invy+ynum*3,-1/3*kbond); % Coefficient for vy(i,j+1)
                    R(invy,1)=0;
                elseif(bc==NEUMANN_BC || bc == PSI_BC)
                    % Neumann dvy/dx=0: vy(i,j)-vy(i,j+1)=0
                    %                 eqv. for psi boundary condtion
                    Lf(invy,invy,1*kbond);                   % Coefficient for vy(i,j)
                    Lf(invy,invy+ynum*3,-1*kbond);                   % Coefficient for vy(i,j+1)
                    R(invy,1)=0;
                end
            end
            % Right boundary, inner points (j=xnum-1, 1<i<ynum)
            if(j==xnum-1 && i>1 && i<ynum)
                if(bc==DIRICHLET_BC)
                    %                 Dirichlet or no slip: vy=0: vy(i,j)-1/3*vy(i,j-1)=0
                    Lf(invy,invy,1*kbond); % Coefficient for vy(i,j)
                    Lf(invy,invy-ynum*3,-1/3*kbond); % Coefficient for vy(i,j-1)
                    R(invy,1)=0;
                elseif(bc==NEUMANN_BC || bc==PSI_BC)
                    % Free slip dvy/dx=0: vy(i,j)-vy(i,j-1)=0
                    %                 eqv. to psi boundary condition
                    Lf(invy,invy,1*kbond);                    % Coefficient for vy(i,j)
                    Lf(invy,invy-ynum*3,-1*kbond);                   % Coefficient for vy(i,j-1)
                    R(invy,1)=0;
                end
            end
            %Internal nodes: mu(d2vy/dx2+d2vy/dy2)-dP/dy=-gy*RHO
        else
            %dSyy/dy=2*muc(i+1,j+1)*(vy(i+1,j)-vy(i,j))/dy^2-2*muc(i,j+1)*(vy(i,j)-vy(i-1,j))/dy^2
            Lf(invy, invy+3, 2*muc(i+1,j+1)/ystp^2);                 % Coefficient for vy(i+1,j)
            Lf(invy, invy-3, 2*muc(i,j+1)/ystp^2);                   % Coefficient for vy(i-1,j)
            Lf(invy, invy, -2*muc(i+1,j+1)/ystp^2-2*muc(i,j+1)/ystp^2); % Coefficient for vy(i,j)
            
            %dSxy/dx=mus(i,j+1)*((vy(i,j+1)-vy(i,j))/dx^2+(vx(i,j+1)-vx(i-1,j+1))/dx/dy)-
            %         -mus(i,j)*((vy(i,j)-vy(i,j-1))/dx^2+(vx(i,j)-vx(i-1,j))/dx/dy)-
            Lf(invy,invy+ynum*3, mug(i,j+1)/xstp^2);                     % Coefficient for vy(i,j+1)
            Lf(invy,invy-ynum*3, mug(i,j)/xstp^2);                       % Coefficient for vy(i,j-1)
            Lf(invy,invy, -mug(i,j+1)/xstp^2-mug(i,j)/xstp^2); % ADD coefficient for vy(i,j)
            Lf(invy,invx+ynum*3, mug(i,j+1)/xstp/ystp);                  % Coefficient for vx(i,j+1)
            Lf(invy,invx+ynum*3-3, -mug(i,j+1)/xstp/ystp);                 % Coefficient for vx(i-1,j+1)
            Lf(invy,invx, -mug(i,j)/xstp/ystp);                   % Coefficient for vx(i,j)
            Lf(invy,invx-3, mug(i,j)/xstp/ystp);                    % Coefficient for vx(i-1,j)
            
            % -dP/dy=(P(i,j+1)-P(i+1,j+1))/dx
            Lf(invy,inp+ynum*3, kcont/ystp);                             % Coefficient for P(i,j+1)
            Lf(invy,inp+3+ynum*3, -kcont/ystp);                            % Coefficient for P(i+1,j+1)
            % Right part: da/dy = (a(i+1,j+1)-a(i,j+1))/dy
            R(invy,1)               = ((muc(i+1,j+1)+muc(i,j+1))/2 + (lambdac(i+1,j+1)+lambdac(i,j+1))/2)*(a(i+1,j+1)-a(i,j+1))/ystp;
        end
        
    end
end



% Constraints when neumann boundary for displacement is used for all the boundaries
if (bc==NEUMANN_BC)
    % For rigid body translation
    % % Let's set the average of the each of the displacement components to be
    % % zero!
    for i = 2:ynum-1
        for j = 2:xnum-1
            invx = ((j-1)*ynum+i)*3 - 1;
            invy = invx + 1;
            Lf(1+xnum*ynum*3,invx, 1*kbond); % Coefficient for vx(i+1,j)
            Lf(2+xnum*ynum*3,invy, 1*kbond); % Coefficient for vx(i-1,j)
        end
    end
    R(1+xnum*ynum*3,1) = 0;
    R(2+xnum*ynum*3,1) = 0;
    
    
    % For rigid body rotation
    % At all interior points, mean of rigid body rotation tensor components set
    % to zero:
    % i.e. dvx/dy - dvy/dx = 0.
    % dvx/dy: (vx(i+1,j)-vx(i-1,j))/(2*ystp) - (vy(i,j+1)-vy(i,j-1))/(2*ystp) = 0;
    for i = 2:ynum-1
        for j = 2:xnum-1
            invx = ((j-1)*ynum+i)*3 - 1;
            invy = invx + 1;
            Lf(3+xnum*ynum*3,invx + 3, kbond/(2*xstp)); % Coefficient for vx(i+1,j)
            Lf(3+xnum*ynum*3,invx - 3, -kbond/(2*xstp)); % Coefficient for vx(i-1,j)
            
            Lf(3+xnum*ynum*3,invy + ynum*3, -kbond/ystp); % Coefficient for vy(i,j+1)
            Lf(3+xnum*ynum*3,invy - ynum*3, kbond/ystp); % Coefficient for vy(i,j-1)
        end
    end
    R(3+xnum*ynum*3,1) = 0;
    ecnstrnt_indx = 4;
    
    % If single point is considered instead of for all interior points:
    % Fix the center, i.e. the vx and vy corresponding to (ynum/2,xnum/2):
    % i.e. vx(ynum/2,xnum/2) = 0;
    % and vy(ynum/2,xnum/2) = 0;
    
    
    
    
    % top left corner
    % i = 2;
    % j = 2;
    
    % top middle
    % i = 2;
    % j = ceil(xnum/2);
    
    % top right corner
    % i = 2;
    % j = xnum-1;
    
    % center
    % i = ceil(ynum/2);
    % j = ceil(xnum/2);
    
    % bottom left corner
    % i = ynum - 1;
    % j = 2;
    
    % bottom middle
    % i = ynum - 1;
    % j = ceil(xnum/2);
    
    % bottom right
    % i = ynum - 1;
    % j = xnum - 1;
    
    % invx = ((j-1)*ynum+i)*3-1;
    % invy = invx + 1;
    %
    % % % vx(i,j) = 0;
    % L(1+xnum*ynum*3,invx) = 1*kbond;
    % R(1+xnum*ynum*3,1) = 0;
    % %
    % % % vy(i,j) = 0;
    % L(2+xnum*ynum*3,invy) = 1*kbond;
    % R(2+xnum*ynum*3,1) = 0;
    
    
    
    % Now set dvx/dy - dvy/dx = 0 as another constraint!
    % dvx/dy: (vx(i+1,j)-vx(i-1,j))/(2*ystp) - (vy(i,j+1)-vy(i,j-1))/(2*ystp) = 0;
    % At single point
    % L(3+xnum*ynum*3,invx + 3       )    =    kcont/(2*xstp); % Coefficient for vx(i+1,j)
    % L(3+xnum*ynum*3,invx - 3)    =   -kcont/(2*xstp); % Coefficient for vx(i-1,j)
    % % % % %
    % % % % %-dvy/dx: -(vy(i,j+1)+vy(i,j-1))/(2*xstp)
    % L(3+xnum*ynum*3,invy + ynum*3  )    =   -kcont/ystp; % Coefficient for vy(i,j+1)
    % L(3+xnum*ynum*3,invy - ynum*3)    =   kcont/ystp; % Coefficient for vy(i,j-1)
    %
    % R(3+xnum*ynum*3,1) = 0;
    
    if (pChoice == ENFORCE_CC) %Since CC is enforced, extra constraint needed here to fix pressure value to make it full rank!
        inp     =   ((p0cellx-1)*ynum+p0celly)*3-2;
        Lf(ecnstrnt_indx + xnum*ynum*3,inp, 1*kbond);  %Coeff. for p0cell
        R(ecnstrnt_indx + xnum*ynum*3,1) = p0cell;
        ecnstrnt_indx = ecnstrnt_indx + 1;
    end
end

if (bc ~= NEUMANN_BC && pChoice == ENFORCE_CC)
    inp     =   ((p0cellx-1)*ynum+p0celly)*3-2;
    Lf(1+xnum*ynum*3,inp, 1*kbond);  %Coeff. for p0cell
    R(1+xnum*ynum*3,1) = p0cell;
    ecnstrnt_indx = ecnstrnt_indx + 1;
end


% Create sparse matrix from the indices and values
L = sparse(i_indx(1:vals_cntr),j_indx(1:vals_cntr),L_vals(1:vals_cntr),ecnstrnt+xnum*ynum*3,xnum*ynum*3);

toc;
display('matrix built');

% write the matrix to solve with PetsC:
% writeAandBtoFile(L,R,'petsC/A.txt','petsC/b.txt');
tic;

%Obtaining vector of solutions S()
% S=L\R;

tol = 1e-12; maxit = 40;
[LL,UU] = ilu(L,struct('type','ilutp','droptol',1e-10));
% [S,fl1,rr1,it1,rv1] = bicgstab(L,R,tol,maxit,LL,UU);
[S,fl1,rr1,it1,rv1] = gmres(L,R,100,tol,maxit,LL,UU);


% options = AMGinit(L);
% [prec options] = AMGfactor(L,options);
% [S, options] = AMGsolver(S,prec,options,R);


toc;
display('system solver timed');
% Reload solutions to 2D p(), vx(), vy() arrays
% Dimensions of arrays are reduced compared to the basic grid
p=zeros(ynum,xnum);
vy=zeros(ynum,xnum);
vx=zeros(ynum,xnum);
% Process all Grid points
for i=1:1:ynum
    for j=1:1:xnum
        % Global index for P, vx, vy in S()
        inp=((j-1)*ynum+i)*3-2; % P
        invx=inp+1;
        invy=inp+2;
        % P
        p(i,j)=S(inp)*kcont;
        % vx
        vx(i,j)=S(invx);
        % vy
        vy(i,j)=S(invy);
    end
end

% save('results2D.mat','vx','vy','p','L','R');

% Compute vx,vy for internal nodes
% Note here too, vx1 and vy1 will have dummy 1st row and 1st col with zero!
% So compatible with the dummy atrophy!
vx1=zeros(ynum,xnum);
vy1=zeros(ynum,xnum);
% Process internal Grid points
for i=2:1:ynum-1
    for j=2:1:xnum-1
        % vx
        vx1(i,j)=(vx(i-1,j)+vx(i,j))/2;
        % vy
        vy1(i,j)=(vy(i,j-1)+vy(i,j))/2;
    end
end



%Plotting solution
% Making new figure
figure;

% Plotting pressure as colormap
% subplot(1,2,1);
% pcolor(xc/1000,yc/1000,p(2:1:ynum,2:1:xnum));%*1e-9);      % making a colormap
% shading interp;     % making smooth transitions between colors
imagesc(p(2:1:ynum,2:1:xnum));%*1e-9);
colorbar;           % showing a colorbar for the map
hold on;            % continuing plotting on the colormap
% Plotting velocity vector as arrows using internal nodes only
% quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
% quiver(vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
quiver(vx1(2:1:ynum-1,2:1:xnum-1)*10000,vy1(2:1:ynum-1,2:1:xnum-1)*10000,'k','Autoscale','Off'); % making field of arrows
hold off;           % stop plotting on the colormap
box on;             % making a box around the plot
title('Pressure (color), velocity (arrows)'); % title for the plot
% xlabel('x, m');        % title for the horizontal axis
% ylabel('y, m');        % title for the vertical axis
axis ij image ;     % directing vertical axis downward, making proper dimensions
% axis([0 xsize/1000 0 ysize/1000]); % Making axes limits


% Plotting viscosity
figure,
imagesc(mug), title(['muBrain: ' num2str(muBrain) ' muCSF: ' num2str(muCSF)]);
axis equal;
% Plotting density as colormap
% subplot(1,2,2);
% pcolor(x/1000,y/1000,rho);      % making a colormap
% shading interp;     % making smooth transitions between colors
% colorbar;           % showing a colorbar for the map
% hold on;            % continuing plotting on the colormap
% Plotting velocity vector as arrows using internal nodes only
% quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
% hold off;           % stop plotting on the colormap
% box on;             % making a box around the plot
% title('Density (color, kg/m^3), velocity (arrows)');   % title for the plot
% xlabel('x, km');        % title for the horizontal axis
% ylabel('y, km');        % title for the vertical axis
% axis ij image ;     % directing vertical axis downward, making proper dimensions
% axis([0 xsize/1000 0 ysize/1000]); % Making axes limits

figure,
subplot(121),
% Display atrophy except the dummy ones (i.e. 1st row and 1st column)
imagesc(a(2:end,2:end)), title('atrophy, a');
% axis image;
axis equal;

soltn(:,:,1) = vx;
soltn(:,:,2) = vy;
div_ans = div2D_new(soltn,xstp,ystp);
% Last row/col contains false div. so remove it
div_ans = div_ans(1:end-1,1:end-1);
% Now add dummy on first row and col to make it comparible with a:
div_ans = [zeros(1,xnum); zeros(ynum-1,1) div_ans];
% Display divergence of the solution except for the dummy row and col.
subplot(122), imagesc(div_ans(2:end,2:end)), title('divergence of the solution displ. field ');
axis equal;


% clims = [min([min(a) min(div_ans)]) max([max(a) max(div_ans)])];
figure,
% subplot(121), imagesc(a,clims), title('atrophy');
% subplot(122), imagesc(-1*(div_ans),clims), title('divergence');
% subplot(121), imagesc(a), title('atrophy');
% subplot(122), imagesc(div_ans), title('divergence');
imagesc(a+div_ans), title('constraint');
axis equal;

% Re-add the removed part with zero velocities.
% Note that in using vx1, vx1 already has a zero first row and col due to
% dummy addition during interpolations from vx before this concatenation.
% Always expects a seg_mask with CSF in the periphery
    [m n] = size(div_slice);
    vc = size(vx1,2);
    vx2 = cat(1,zeros(min_max_row_col(1)-extra_CSF-1-1, vc),vx1,zeros(m-min_max_row_col(2)-extra_CSF-1+1,vc));
    vx2 = cat(2,zeros(m, min_max_row_col(3)-extra_CSF-1-1),vx2,zeros(m, n-min_max_row_col(4)-extra_CSF-1+1));
    
    vy2 = cat(1,zeros(min_max_row_col(1)-extra_CSF-1-1, vc),vy1,zeros(m-min_max_row_col(2)-extra_CSF-1+1,vc));
    vy2 = cat(2,zeros(m, min_max_row_col(3)-extra_CSF-1-1),vy2,zeros(m, n-min_max_row_col(4)-extra_CSF-1+1));
    
    v2(:,:,1) = vx2;
    v2(:,:,2) = vy2;
end