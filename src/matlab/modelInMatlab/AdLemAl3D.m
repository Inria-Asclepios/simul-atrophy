% 3D Solution of the Brain deformation with irregular domain (or brain-CSF
% considered as an interface and extending the boundary to the skull.
% Regular staggered grid using a pressure-velocity formulation
% for a medium with varying lame's parameters.

function AdLemAl3D(muBrain,muRatio,lambdaBrain,lambdaRatio,seg_mask,div_map,USE_REAL_DATA)

if(nargin < 7)
    USE_REAL_DATA = false;
end


if (USE_REAL_DATA == true)
   % Since mu and lambda Brain are greater than CSF, the ratio muB/muC given
    % as input in the form muRatio should be greater or equal to 1.
    if(lambdaRatio < 1 || muRatio < 1)
        error('invalid muRatio or invalid lamdaRatio, they should be greater than one');
    end

%     load('/home/bkhanal/staggered3D/augLagr/real_data.mat');
    [brain_mask min_max_yxz] = getCuboidConvexHull(seg_mask);
    [yn xn zn] = size(brain_mask);
    a = double(div_map(min_max_yxz(1) : min_max_yxz(2),...
        min_max_yxz(3) : min_max_yxz(4), min_max_yxz(5) : min_max_yxz(6)));
    
    %pad with extra_CSF width of zeros all around.
    eCSF = 2;
    brain_mask = padarray(brain_mask,[eCSF eCSF eCSF]);
    
    % distributed source outside convex hull(padded areas), 
    % corners where CC is not enforced! Total pixels of padded areas of
    % width w is (xn+2w)(yn+2w)(zn+2w) - xn.yn.zn
    a = padarray(a,[eCSF eCSF eCSF],-sum(a(:))/...
        ((xn+2*eCSF)*(yn+2*eCSF)*(zn+2*eCSF) - xn*yn*zn));
    
    %Grid nodes have one dimension greater than the cell centers!
    yn = yn+eCSF*2+1;
    xn = xn+eCSF*2+1;
    zn = zn+eCSF*2+1;
    
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
    xn    =   11;             % Horizontal
    yn    =   11;             % Vertical
    zn    =   11;             % Horizontal towards the screen
    
    % Call a function:
    center_mask = [round(xn/2) round(yn/2) round(zn/2)];
    % brain_mask = sphereMask(center_mask,xnum/4,xnum-1,ynum-1);
    % brain_mask = ellipsoidMask(center_mask,xnum/3,ynum/4,znum/4,xnum-1,ynum-1,znum-1);
    % brain_mask = ellipsoidMask(center_mask,7,4,4,xnum-1,ynum-1,znum-1);
    brain_mask = sphereMask(center_mask,floor(xn/2),xn-1,yn-1,zn-1);
    
    % brain_mask = sphereMask(center_mask,7,xnum-1,ynum-1,znum-1);
    
    % Get the atrophy data at the cell-centers
    % atrophy at the center of the cells
    a = computeAtrophy3D(yn-1,xn-1,zn-1);
    
    % distributed source outside brain, excluding 12 edges where CC is not
    % enforced!
%     Before dummy, size of total size of a is (xn-1)*(yn-1)*(zn-1)!!
    a(brain_mask==false) = -sum(a(:))/((xn-1)*(yn-1)*(zn-1) - sum(brain_mask(:)));
    % add a dummy atrophy corresponding to dummy pressure values!
    a = padarray(a,[1 1 1],'pre');
%     a = cat(3,zeros(yn-1,xn-1,1),a);
%     a = [zeros(1,xn,zn); zeros(yn-1,1,zn) a];
    
end

NS     =  xn*yn*zn*3;     % Total number of variables to compute values of.

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

% Now Add dummy (duplicate) on 1st row,col and plane of centered cells:
muc = padarray(muc,[1 1 1],'replicate','pre');

% Now Add dummy duplicate on 1st row,col and plane of centered cells:
lambdac = padarray(lambdac,[1 1 1],'replicate','pre');

% Adaptive Lagrange parameter:
rc = ones(yn-1,xn-1,zn-1);
rc(brain_mask==true) = muBrain/10;
rc(brain_mask==false) = muCSF/0.1;

% Compute at the planes:
ry = (rc(1:end-1,:,:) + rc(2:end,:,:))/2;
ry = padarray(ry,[1 0 0],'replicate','symmetric');
ry = padarray(ry,[0 1 1],'replicate','post');

rx = (rc(:,1:end-1,:) + rc(:,2:end,:))/2;
rx = padarray(rx,[0 1 0],'replicate','symmetric');
rx = padarray(rx,[1 0 1],'replicate','post');

rz = (rc(:,:,1:end-1) + rc(:,:,2:end))/2;
rz = padarray(rz,[0 0 1],'replicate','symmetric');
rz = padarray(rz,[1 1 0],'replicate','post');

% Now extend a dummy on the 1st row and 1st col of centered cells:
rc = padarray(rc,[1 1 1], 'replicate','pre');

% Computing Kcont and Kbond coefficients
kbond   =   4*muBrain/(xstp+ystp+zstp)^2;
% kbond   =   4*muCSF/(xstp+ystp+zstp)^2;

DIRICHLET_BC = 0;
NEUMANN_BC = 1;

% Boundary Condition Options:
bc = DIRICHLET_BC; %dirichlet condition on displacement.
% bc = NEUMANN_BC; %neumann condition on displacement.

% Number of Extra constraints required, 3 for translation, 3 for rotation
% if (bc == NEUMANN_BC)
if (bc == NEUMANN_BC)
    ecnstrnt = 6;
else
    ecnstrnt = 0;
end


% index keeping track of the ecnstrnt
ecnstrnt_indx = 1;


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
% R       =   zeros(ecnstrnt+ NS,1);

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



    function R = setRhs(p_in)
        R = zeros(ecnstrnt + NS,1);
        for ii=1:1:yn
            for jj=1:1:xn
                for kk = 1:1:zn
                    
                    % Global index for P, vx, vy and vz in the current node
                    i_nvx    =   ((kk-1)*xn*yn + (jj-1)*yn+ii)*3-2; % invx
                    i_nvy    =   i_nvx+1;
                    i_nvz    =   i_nvx+2;
                    
                    
                    % x-Stokes equation
                    % Ghost vx unknowns (i=yn || k=zn) and boundary nodes (i=1, i=yn-1, j=1, j=xn, k=1, k=zn-1)
                    if(ii==yn || kk==zn || ii==1 || ii==yn-1 || jj==1 || jj==xn || kk==1 || kk==zn-1)
                        % Ghost vx unknowns (i=yn or k=zn: vx(i,j,k)=0
                        if(ii==yn || kk==zn)
                            R(i_nvx,1   )        =   0;          % Right-hand-side part
                        end
                        
                        % Left and Right boundaries (j=1, j=xn)
                        %             Left boundary
                        if((jj==1 ) && (ii<yn) && (kk<zn))
                            if (bc == DIRICHLET_BC)
                                %                 Dirichlet condition: vx(i,j,k) = 0
                                R(i_nvx,1)   = 0;
                            elseif (bc == NEUMANN_BC)
                                % Neumann condition: -vx(i,j+1,k)+vx(i,j,k)=0
                                R(i_nvx,1)           =   0;          % Right-hand-side part
                                %
                                %                 elseif (bc == PSI_BC)
                                %                     %                 psi boudnary condition
                                %                     L(i_nvx,i_nvx) = 1*kbond;
                                %                     R(i_nvx,1) = kbond*(psi(i,j)/xstp);
                            end
                        end
                        %             Right boundary
                        if((jj==xn) && (ii<yn) && (kk<zn))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition: vx(i,j,k) = 0
                                R(i_nvx,1)    = 0;
                            elseif(bc==NEUMANN_BC)
                                % %                 Neumann condition: vx(i,j,k) - vx(i,j-1,k) = 0
                                R(i_nvx,1) = 0;
                                %                 elseif(bc==PSI_BC)
                                %                     %             psi boundary condition
                                %                     L(i_nvx,i_nvx) = 1*kbond;
                                %                     R(i_nvx,1) = kbond*(psi(i,j)/xstp);
                            end
                        end
                        %             Upper and lower boundaries
                        % Upper boundary, inner points w.r.t x-axis, (i=1, k<zn, 1<j<xn)
                        if(ii==1 && kk<zn && jj>1 && jj<xn)
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition: vx(i,j,k) = 0
                                % % equivalent to No slip vx=0: vx(i,j,k)-1/3*vx(i+1,j,k)=0.
                                %                     this extrapolation assumption might have affect on
                                %                     not matching the divergence theorem exactly!!
                                R(i_nvx,1)          =   0;          % Right-hand-side part
                            elseif(bc==NEUMANN_BC)
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvx/dy=0: vx(i+1,j,k)-vx(i,j,k)=0
                                %                 Same for psi boundary condition!
                                R(i_nvx,1)           =   0;          % Right-hand-side part
                            end
                        end
                        % Lower boundary, inner points w.r.t x-axis (i=yn-1, k<zn, 1<j<xn)
                        if(ii==yn-1 && kk<zn && jj>1 && jj<xn)
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet or No slip:
                                % % No slip vx=0: vx(i,j,k)-1/3*vx(i-1,j,k)=0
                                R(i_nvx,1)         =   0; % Right part
                            elseif(bc==NEUMANN_BC )
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvx/dy=0: vx(i,j)-vx(i-1,j)=0
                                %                 of psi condition
                                R(i_nvx,1)           =   0; % Right-hand-side part
                            end
                        end
                        
                        
                        %             Front and back boundaries
                        % Front boundary, inner points w.r.t, (k=1, 1<j<xn, 1<i<yn-1)
                        if(kk==1 && jj>1 && jj<xn && ii>1 && (ii<(yn-1)))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition: vx(i,j,k) = 0
                                % % equivalent to No slip vx=0: vx(i,j,k)-1/3*vx(i,j,k+1)=0.
                                %                     this extrapolation assumption might have affect on
                                %                     not matching the divergence theorem exactly!!
                                R(i_nvx,1)          =   0;          % Right-hand-side part
                            elseif(bc==NEUMANN_BC)
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvx/dy=0: vx(i,j,k+1)-vx(i,j,k)=0
                                %                 Same for psi boundary condition!
                                R(i_nvx,1)           =   0;          % Right-hand-side part
                            end
                        end
                        % Back boundary, inner points (k=zn-1, 1<j<xnum, 1<i<yn-1)
                        if(kk==zn-1 && jj>1 && jj<xn && ii>1 && (ii<(yn-1)))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet or No slip:
                                % % No slip vx=0: vx(i,j,k)-1/3*vx(i,j,k-1)=0
                                R(i_nvx,1)         =   0; % Right part
                            elseif(bc==NEUMANN_BC )
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvx/dy=0: vx(i,j,k)-vx(i,j,k-1)=0
                                %                 of psi condition
                                R(i_nvx,1)           =   0; % Right-hand-side part
                            end
                        end
                        
                        
                        %Internal nodes: mu(d2vx/dx2+d2vx/dy2)-dP/dx=0
                    else
                        % Right part:dp/dx + (-r+mu+lambda)da/dx
                        R(i_nvx,1) = (p_in(ii+1,jj+1,kk+1) - p_in(ii+1,jj,kk+1))/xstp...
                            + (-rx(ii,jj,kk) + (muc(ii+1,jj+1,kk+1)+muc(ii+1,jj,kk+1))/2 ...
                            + (lambdac(ii+1,jj+1,kk+1)+lambdac(ii+1,jj,kk+1))/2)...
                            * (a(ii+1,jj+1,kk+1)-a(ii+1,jj,kk+1))/xstp;
                    end
                    
                    % y-Stokes equation
                    % Ghost vy unknowns (j=xn || k=zn) and boundary nodes (i=1, i=yn, j=1, j=xn-1, k=1, k=zn-1)
                    if(jj==xn || kk==zn || ii==1 || ii==yn || jj==1 || jj==xn-1 || kk==1 || kk==zn-1)
                        % Ghost vy unknowns (j=xn or k=zn: vy(i,j)=0
                        if(jj==xn || kk==zn)
                            R(i_nvy,1)           =   0;
                        end
                        % Upper boundary, i=1
                        if((ii==1) && (jj<xn) && (kk<zn))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition vy = 0
                                R(i_nvy,1) = 0;
                            elseif(bc==NEUMANN_BC)
                                % %            Neumann codition: -vy(i+1,j,k) + vy(i,j,k) = 0;
                                R(i_nvy,1)           =   0;
                                %                 elseif(bc==PSI_BC)
                                %                     %             psi boundary condtion:
                                %                     L(i_nvy,i_nvy) = 1*kbond;
                                %                     R(i_nvy,1) = kbond*(psi(i,j)/ystp);
                            end
                        end
                        %             Lower boundary, i=yn
                        if(ii==yn && jj<xn && kk<zn)
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition vy = 0
                                R(i_nvy,1) = 0;
                            elseif(bc==NEUMANN_BC)
                                % %           Neumann codition: vy(i,j,k) - vy(i-1,j,k) = 0;
                                R(i_nvy,1) = 0;               %Right side
                                %                 elseif(bc==PSI_BC)
                                %                     %                 psi boundary condtion:
                                %                     L(i_nvy,i_nvy) = 1*kbond;
                                %                     R(i_nvy,1) = kbond*(psi(i,j)/ystp);
                            end
                        end
                        
                        %             Left and right boundaries
                        % Left boundary, inner points w.r.t. y-axis, (j=1, k<zn, 1<i<yn,)
                        if(jj==1 && kk<zn && ii>1 && ii<yn)
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet or eqv. to no slip:
                                %             % No slip vy=0: vy(i,j,k)-1/3*vy(i,j+1,k)=0
                                R(i_nvy,1)=0;
                            elseif(bc==NEUMANN_BC)
                                %                     elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann dvy/dx=0: vy(i,j,k)-vy(i,j+1,k)=0
                                %                 eqv. for psi boundary condtion
                                R(i_nvy,1)=0;
                            end
                        end
                        % Right boundary, inner points w.r.t y-axis, (j=xn-1, k<zn, 1<i<yn)
                        if(jj==xn-1 && kk<zn && ii>1 && ii<yn)
                            if(bc==DIRICHLET_BC)
                                %                 Dirichlet or no slip: vy=0: vy(i,j,k)-1/3*vy(i,j-1,k)=0
                                R(i_nvy,1)=0;
                            elseif(bc==NEUMANN_BC)
                                %                 elseif(bc==NEUMANN_BC || bc==PSI_BC)
                                % Free slip dvy/dx=0: vy(i,j,k)-vy(i,j-1,k)=0
                                %                 eqv. to psi boundary condition
                                R(i_nvy,1)=0;
                            end
                        end
                        
                        %             Front and back boundaries
                        % Front boundary, inner points (k=1, 1<j<xn-1, 1<i<yn)
                        if(kk==1 && ii>1 && ii<yn && jj>1 && (jj<(xn-1)))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition: vy(i,j,k) = 0
                                % % equivalent to No slip vy=0: vy(i,j,k)-1/3*vy(i,j,k+1)=0.
                                %                     this extrapolation assumption might have affect on
                                %                     not matching the divergence theorem exactly!!
                                R(i_nvy,1)          =   0;          % Right-hand-side part
                            elseif(bc==NEUMANN_BC)
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvy/dz=0: vy(i,j,k)-vy(i,j,k+1)=0
                                %                 Same for psi boundary condition!
                                R(i_nvy,1)           =   0;          % Right-hand-side part
                            end
                        end
                        % Back boundary, inner points (k=zn-1, 1<j<xn-1, 1<i<yn)
                        if(kk==zn-1 && ii>1 && ii<yn && jj>1 && (jj<(xn-1)))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet or No slip:
                                % % No slip vy=0: vy(i,j,k)-1/3*vy(i,j,k-1)=0
                                R(i_nvy,1)         =   0; % Right part
                            elseif(bc==NEUMANN_BC )
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvy/dz=0: vy(i,j,k)-vy(i,j,k-1)=0
                                %                 of psi condition
                                R(i_nvy,1)           =   0; % Right-hand-side part
                            end
                        end
                        
                        
                        %Internal nodes:
                        %mu(d2vy/dx2+d2vy/dy2+d2vy/dz2)-dP/dy=(mu+lambda)da/dy
                    else
                        % dp/dy + (-r + mu+lambda)da/dy
                        R(i_nvy,1) = (p_in(ii+1,jj+1,kk+1) - p_in(ii,jj+1,kk+1))/xstp + ...
                            + (-ry(ii,jj,kk)+(muc(ii+1,jj+1,kk+1)+muc(ii,jj+1,kk+1))/2 ...
                            + (lambdac(ii+1,jj+1,kk+1)+lambdac(ii,jj+1,kk+1))/2)...
                            *(a(ii+1,jj+1,kk+1)-a(ii,jj+1,kk+1))/ystp;
                    end
                    
                    
                    % z-Stokes equation
                    % Ghost vz unknowns (i=yn || j=xn) and boundary nodes (i=1, i=yn-1, j=1, j=xn-1, k=1, k=zn)
                    if(ii==yn || jj==xn ||  kk==1 || kk==zn || jj==1 || jj==xn-1 || ii==1 || ii==yn-1)
                        % Ghost vz unknowns (j=xn or i=yn: vz(i,j,k)=0
                        if( ii==yn || jj==xn)
                            R(i_nvz,1)           =   0;
                        end
                        % Front boundary
                        if((kk==1) && (ii<yn) && (jj<xn))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition vz = 0
                                R(i_nvz,1) = 0;
                            elseif(bc==NEUMANN_BC)
                                % %            Neumann codition: -vz(i,j,k+1) + vy(i,j,k) = 0;
                                R(i_nvz,1)           =   0;
                                %                 elseif(bc==PSI_BC)
                                %                     %             psi boundary condtion:
                                %                     L(i_nvy,i_nvy) = 1*kbond;
                                %                     R(i_nvy,1) = kbond*(psi(i,j)/ystp);
                            end
                        end
                        %  Back boundary, k=zn
                        if(kk==zn && ii<yn && jj<xn)
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition vz = 0
                                R(i_nvz,1) = 0;
                            elseif(bc==NEUMANN_BC)
                                % %           Neumann codition: vz(i,j,k) - vz(i,j,k-1) = 0;
                                R(i_nvz,1) = 0;               %Right side
                            end
                        end
                        
                        %             Left and right boundaries
                        % Left boundary, inner points w.r.t. z-axis, (j=1, i<yn, 1<k<zn,)
                        if(jj==1 && ii<yn && kk>1 && kk<zn)
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet or eqv. to no slip:
                                %             % No slip vz=0: vz(i,j,k)-1/3*vz(i,j+1,k)=0
                                R(i_nvz,1)=0;
                            elseif(bc==NEUMANN_BC)
                                %                     elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann dvz/dx=0: vz(i,j,k)-vz(i,j+1,k)=0
                                %                 eqv. for psi boundary condtion
                                R(i_nvz,1)=0;
                            end
                        end
                        % Right boundary, inner points w.r.t z-axis, (j=xn-1, i<yn, 1<k<zn)
                        if(jj==xn-1 && ii<yn && kk>1 && kk<zn)
                            if(bc==DIRICHLET_BC)
                                %                 Dirichlet or no slip: vz=0: vz(i,j,k)-1/3*vz(i,j-1,k)=0
                                R(i_nvz,1)=0;
                            elseif(bc==NEUMANN_BC)
                                %                 elseif(bc==NEUMANN_BC || bc==PSI_BC)
                                % Free slip dvz/dx=0: vz(i,j,k)-vz(i,j-1,k)=0
                                %                 eqv. to psi boundary condition
                                R(i_nvz,1)=0;
                            end
                        end
                        
                        %   Upper and Lower Boundaries:
                        % Upper boundary, inner points (i=1, 1<j<xn-1, 1<k<zn)
                        if(ii==1 && kk>1 && kk<zn && jj>1 && (jj<(xn-1)))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet condition: vz(i,j,k) = 0
                                % % equivalent to No slip vz=0: vz(i,j,k)-1/3*vz(i+1,j,k)=0.
                                %                     this extrapolation assumption might have affect on
                                %                     not matching the divergence theorem exactly!!
                                R(i_nvz,1)          =   0;          % Right-hand-side part
                            elseif(bc==NEUMANN_BC)
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvz/dy=0: vz(i+1,j,k)-vz(i,j,k)=0
                                %                 Same for psi boundary condition!
                                R(i_nvz,1)           =   0;          % Right-hand-side part
                            end
                        end
                        % Lower boundary, inner points (i=yn-1, 1<j<xn-1, 1<k<zn)
                        if((ii==(yn-1)) && kk>1 && kk<zn && jj>1 && (jj<(xn-1)))
                            if(bc==DIRICHLET_BC)
                                %                     Dirichlet or No slip:
                                % % No slip vz=0: vz(i,j,k)-1/3*vz(i-1,j,k)=0
                                R(i_nvz,1)         =   0; % Right part
                            elseif(bc==NEUMANN_BC )
                                %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                                % Neumann condition dvz/dy=0: vz(i,j,k)-vz(i-1,j,k)=0
                                %                 of psi condition
                                R(i_nvz,1)           =   0; % Right-hand-side part
                            end
                        end
                        
                        %Internal nodes:
                        %mu(d2vz/dx2+d2vz/dy2+d2vz/dz2)-dP/dz=(mu+lambda)da/dz
                    else

                        % Right part: dp/dz + (-r+mu+lambda)*da/dz
                        % 
                        R(i_nvz,1) = (p_in(ii+1,jj+1,kk+1)-p_in(ii+1,jj+1,kk))/zstp...
                            + (-rz(ii,jj,kk)+(muc(ii+1,jj+1,kk+1)+muc(ii+1,jj+1,kk))/2 ...
                            + (lambdac(ii+1,jj+1,kk+1)+lambdac(ii+1,jj+1,kk))/2)...
                            *(a(ii+1,jj+1,kk+1)-a(ii+1,jj+1,kk))/zstp;
                    end
                end
            end
        end
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
            invx    =   ((k-1)*xn*yn + (j-1)*yn+i)*3-2; % invx
            invy    =   invx+1;
            invz    =   invx+2;
            
            % x-Stokes equation
            % Ghost vx unknowns (i=yn || k=zn) and boundary nodes (i=1, i=yn-1, j=1, j=xn, k=1, k=zn-1)
            if(i==yn || k==zn || i==1 || i==yn-1 || j==1 || j==xn || k==1 || k==zn-1)
                % Ghost vx unknowns (i=yn or k=zn: vx(i,j,k)=0
                if(i==yn || k==zn)
                    Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
%                     R(invx,1   )        =   0;          % Right-hand-side part
                end
                
                % Left and Right boundaries (j=1, j=xn)
                %             Left boundary
                if((j==1 ) && (i<yn) && (k<zn))
                    if (bc == DIRICHLET_BC)
                        %                 Dirichlet condition: vx(i,j,k) = 0
                        Lf(invx,invx, 1*kbond);  %Coeff. for vx(i,j,k)
%                         R(invx,1)   = 0;
                    elseif (bc == NEUMANN_BC)
                        % Neumann condition: -vx(i,j+1,k)+vx(i,j,k)=0
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx+yn*3, -1*kbond);    % coeff. for vx(i,j+1,k)
%                         R(invx,1)           =   0;          % Right-hand-side part
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
%                         R(invx,1)    = 0;
                    elseif(bc==NEUMANN_BC)
                        % %                 Neumann condition: vx(i,j,k) - vx(i,j-1,k) = 0
                        Lf(invx,invx, 1*kbond);  %Coeff. for vx(i,j,k)
                        Lf(invx,invx-yn*3, -1*kbond); %Coeff. for vx(i,j-1,k)
%                         R(invx,1) = 0;
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
                        Lf(invx,invx+3, -1/3*kbond); % Coefficient for vx(i+1,j,k)
%                         R(invx,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i+1,j,k)-vx(i,j,k)=0
                        %                 Same for psi boundary condition!
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx+3, -1*kbond);    % Coefficient for vx(i+1,j,k)
%                         R(invx,1)           =   0;          % Right-hand-side part
                    end
                end
                % Lower boundary, inner points w.r.t x-axis (i=yn-1, k<zn, 1<j<xn)
                if(i==yn-1 && k<zn && j>1 && j<xn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vx=0: vx(i,j,k)-1/3*vx(i-1,j,k)=0
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx-3, -1/3*kbond); % Coefficient for vx(i-1,j,k)
%                         R(invx,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i,j)-vx(i-1,j)=0
                        %                 of psi condition
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx-3, -1*kbond); % Coefficient for vx(i-1,j,k)
%                         R(invx,1)           =   0; % Right-hand-side part
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
                        Lf(invx,invx+ xn*yn*3, -1/3*kbond); % Coefficient for vx(i,j,k+1)
%                         R(invx,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i,j,k+1)-vx(i,j,k)=0
                        %                 Same for psi boundary condition!
                        Lf(invx,invx, 1*kbond);    % Coefficient for vx(i,j,k)
                        Lf(invx,invx + xn*yn*3, -1*kbond);    % Coefficient for vx(i,j,k+1)
%                         R(invx,1)           =   0;          % Right-hand-side part
                    end
                end
                % Back boundary, inner points (k=zn-1, 1<j<xnum, 1<i<yn-1)
                if(k==zn-1 && j>1 && j<xn && i>1 && (i<(yn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vx=0: vx(i,j,k)-1/3*vx(i,j,k-1)=0
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx - xn*yn*3, -1/3*kbond); % Coefficient for vx(i,j,k-1)
%                         R(invx,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvx/dy=0: vx(i,j,k)-vx(i,j,k-1)=0
                        %                 of psi condition
                        Lf(invx,invx, 1*kbond); % Coefficient for vx(i,j,k)
                        Lf(invx,invx - xn*yn*3, -1*kbond); % Coefficient for vx(i,j,k-1)
%                         R(invx,1)           =   0; % Right-hand-side part
                    end
                end
                
                
                %Internal nodes: mu(d2vx/dx2+d2vx/dy2)-dP/dx=0
            else
                %vx coeff.
                Lf(invx, invx+yn*3, (2*muc(i+1,j+1,k+1)/xstp^2)... % Coefficient for vx(i,j+1,k)
                    + rc(i+1,j+1,k+1)/(xstp^2));                              %From d/dx(div(v))
                Lf(invx,invx - yn*3, (2*muc(i+1,j,k+1)/xstp^2)...   
                 + rc(i+1,j,k+1)/(xstp^2));                                 %Coefficient for vx(i,j-1,k)
                Lf(invx,invx, (-2*muc(i+1,j+1,k+1)/xstp^2 - 2*muc(i+1,j,k+1)/xstp^2)...
                    + (-muxy(i+1,j,k)/ystp^2-muxy(i,j,k)/ystp^2)...
                    + (-muxz(i,j,k+1)/zstp^2-muxz(i,j,k)/zstp^2)...
                    -(rc(i+1,j+1,k+1)+rc(i+1,j,k+1))/(xstp^2))                            % Coefficient for vx(i,j,k)
                
                Lf(invx,invx+3, muxy(i+1,j,k)/ystp^2);                             % Coefficient for vx(i+1,j,k)
                Lf(invx,invx-3, muxy(i,j,k)/ystp^2);                               % Coefficient for vx(i-1,j,k)
               
                Lf(invx,invx + xn*yn*3, muxz(i,j,k+1)/zstp^2);                             % Coefficient for vx(i,j,k+1)
                Lf(invx,invx - xn*yn*3, muxz(i,j,k)/zstp^2);                               % Coefficient for vx(i,j,k-1)
                
                %vy coeff.
                
                Lf(invx,invy+3, (muxy(i+1,j,k)/xstp/ystp)...
                    +(rc(i+1,j+1,k+1)/(xstp*ystp)));                                       % Coefficient for vy(i+1,j,k)
                Lf(invx,invy+3-yn*3, (-muxy(i+1,j,k)/xstp/ystp)...
                    -rc(i+1,j,k+1)/(xstp*ystp));                         % Coefficient for vy(i+1,j-1,k)
                Lf(invx,invy, (-muxy(i,j,k)/xstp/ystp)...
                    -rc(i+1,j+1,k+1)/(xstp*ystp));                           % Coefficient for vy(i,j,k)
                Lf(invx,invy-yn*3, (muxy(i,j,k)/xstp/ystp)...
                    +rc(i+1,j,k+1)/(xstp*ystp));                            % Coefficient for vy(i,j-1,k)
                
%                 vz coeff.
                Lf(invx,invz+ xn*yn*3, (muxz(i,j,k+1)/xstp/zstp)...
                    +rc(i+1,j+1,k+1)/(xstp*zstp));                          % Coefficient for vz(i,j,k+1)
                Lf(invx,invz - yn*3 + xn*yn*3, (-muxz(i+1,j,k)/xstp/zstp)...
                    -rc(i+1,j,k+1)/(xstp*zstp));                         % Coefficient for vz(i,j-1,k+1)
                Lf(invx,invz, (-muxz(i,j,k)/xstp/zstp)...
                    -rc(i+1,j+1,k+1)/(xstp*zstp));                           % Coefficient for vz(i,j,k)
                Lf(invx,invz - yn*3, (muxz(i,j,k)/xstp/zstp)...
                    +rc(i+1,j,k+1)/(xstp*zstp));                            % Coefficient for vz(i,j-1,k)

            end
            
            % y-Stokes equation
            % Ghost vy unknowns (j=xn || k=zn) and boundary nodes (i=1, i=yn, j=1, j=xn-1, k=1, k=zn-1)
            if(j==xn || k==zn || i==1 || i==yn || j==1 || j==xn-1 || k==1 || k==zn-1)
                % Ghost vy unknowns (j=xn or k=zn: vy(i,j)=0
                if(j==xn || k==zn)
                    Lf(invy,invy, 1*kbond);                    % Coefficient for vy(i,j,k)
%                     R(invy,1)           =   0;
                end
                % Upper boundary, i=1
                if((i==1) && (j<xn) && (k<zn))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition vy = 0
                        Lf(invy,invy, 1*kbond);
%                         R(invy,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %            Neumann codition: -vy(i+1,j,k) + vy(i,j,k) = 0;
                        Lf(invy,invy, 1*kbond);                    % Coefficient for vy(i,j,k)
                        Lf(invy,invy+3, -1*kbond);    %Coefficient for vy(i+1,j,k)
%                         R(invy,1)           =   0;
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
%                         R(invy,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %           Neumann codition: vy(i,j,k) - vy(i-1,j,k) = 0;
                        Lf(invy,invy, 1*kbond);   %Coeff. for vy(i,j,k)
                        Lf(invy, invy - 3, -1*kbond); %Coeff. for vy(i-1,j,k)
%                         R(invy,1) = 0;               %Right side
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
                        Lf(invy,invy+yn*3, -1/3*kbond); % Coefficient for vy(i,j+1,k)
%                         R(invy,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                     elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann dvy/dx=0: vy(i,j,k)-vy(i,j+1,k)=0
                        %                 eqv. for psi boundary condtion
                        Lf(invy,invy, 1*kbond);                   % Coefficient for vy(i,j,k)
                        Lf(invy,invy+yn*3, -1*kbond);                   % Coefficient for vy(i,j+1,k)
%                         R(invy,1)=0;
                    end
                end
                % Right boundary, inner points w.r.t y-axis, (j=xn-1, k<zn, 1<i<yn)
                if(j==xn-1 && k<zn && i>1 && i<yn)
                    if(bc==DIRICHLET_BC)
                        %                 Dirichlet or no slip: vy=0: vy(i,j,k)-1/3*vy(i,j-1,k)=0
                        Lf(invy,invy, 1*kbond); % Coefficient for vy(i,j,k)
                        Lf(invy,invy-yn*3, -1/3*kbond); % Coefficient for vy(i,j-1,k)
%                         R(invy,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc==PSI_BC)
                        % Free slip dvy/dx=0: vy(i,j,k)-vy(i,j-1,k)=0
                        %                 eqv. to psi boundary condition
                        Lf(invy,invy, 1*kbond);                    % Coefficient for vy(i,j)
                        Lf(invy,invy-yn*3, -1*kbond);                   % Coefficient for vy(i,j-1)
%                         R(invy,1)=0;
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
                        Lf(invy,invy + xn*yn*3, -1/3*kbond); % Coefficient for vy(i,j,k+1)
%                         R(invy,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvy/dz=0: vy(i,j,k)-vy(i,j,k+1)=0
                        %                 Same for psi boundary condition!
                        Lf(invy,invy, 1*kbond);    % Coefficient for vy(i,j,k)
                        Lf(invy,invy + xn*yn*3, -1*kbond);    % Coefficient for vy(i,j,k+1)
%                         R(invy,1)           =   0;          % Right-hand-side part
                    end
                end
                % Back boundary, inner points (k=zn-1, 1<j<xn-1, 1<i<yn)
                if(k==zn-1 && i>1 && i<yn && j>1 && (j<(xn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vy=0: vy(i,j,k)-1/3*vy(i,j,k-1)=0
                        Lf(invy,invy, 1*kbond); % Coefficient for vy(i,j,k)
                        Lf(invy,invy - xn*yn*3, -1/3*kbond); % Coefficient for vy(i,j,k-1)
%                         R(invy,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvy/dz=0: vy(i,j,k)-vy(i,j,k-1)=0
                        %                 of psi condition
                        Lf(invy,invy, 1*kbond); % Coefficient for vy(i,j,k)
                        Lf(invy,invy - xn*yn*3, -1*kbond); % Coefficient for vy(i,j,k-1)
%                         R(invy,1)           =   0; % Right-hand-side part
                    end
                end
                
                
                %Internal nodes:
                %mu(d2vy/dx2+d2vy/dy2+d2vy/dz2)-dP/dy=(mu+lambda)da/dy
            else
                Lf(invy,invy+3, (2*muc(i+1,j+1,k+1)/ystp^2)...
                    +rc(i+1,j+1,k+1)/(ystp*ystp));                 % Coefficient for vy(i+1,j,k)
                Lf(invy,invy-3, (2*muc(i,j+1,k+1)/ystp^2)...
                    +rc(i,j+1,k+1)/(ystp*ystp));                   % Coefficient for vy(i-1,j,k)
                Lf(invy,invy, (-2*muc(i+1,j+1,k+1)/ystp^2-2*muc(i,j+1,k+1)/ystp^2)...
                    +(-muyx(i,j+1,k)/xstp^2-muyx(i,j,k)/xstp^2)...
                    +(-muyz(i,j,k+1)/zstp^2-muyz(i,j,k)/zstp^2)...
                    -rc(i+1,j+1,k+1)/(ystp*ystp) -rc(i,j+1,k+1)/(ystp*ystp)); % Coefficient for vy(i,j,k)
                
                Lf(invy,invy+yn*3, muyx(i,j+1,k)/xstp^2);                     % Coefficient for vy(i,j+1,k)
                Lf(invy,invy-yn*3, muyx(i,j,k)/xstp^2);                       % Coefficient for vy(i,j-1,k)
                
                Lf(invy,invy+xn*yn*3, muyz(i,j,k+1)/zstp^2);                     % Coefficient for vy(i,j,k+1)
                Lf(invy,invy-xn*yn*3, muyz(i,j,k)/zstp^2);                       % Coefficient for vy(i,j,k-1)


                Lf(invy,invx+yn*3, (muyx(i,j+1,k)/xstp/ystp)...
                    +rc(i+1,j+1,k+1)/(ystp*xstp));                  % Coefficient for vx(i,j+1,k)
                Lf(invy,invx-3+yn*3, (-muyx(i,j+1,k)/xstp/ystp)...
                    -rc(i,j+1,k+1)/(ystp*xstp));                 % Coefficient for vx(i-1,j+1,k)
                Lf(invy,invx, (-muyx(i,j,k)/xstp/ystp)...
                    -rc(i+1,j+1,k+1)/(ystp*xstp));                   % Coefficient for vx(i,j,k)
                Lf(invy,invx-3, (muyx(i,j,k)/xstp/ystp)...
                    +rc(i,j+1,k+1)/(ystp*xstp));                    % Coefficient for vx(i-1,j,k)
                
                
                Lf(invy,invz+xn*yn*3, (muyz(i,j,k+1)/zstp/ystp)...
                    +rc(i+1,j+1,k+1)/(ystp*zstp));                  % Coefficient for vz(i,j,k+1)
                Lf(invy,invz-3+xn*yn*3, (-muyz(i,j,k+1)/zstp/ystp)...
                    -rc(i,j+1,k+1)/(ystp*zstp));                 % Coefficient for vz(i-1,j,k+1)
                Lf(invy,invz, (-muyz(i,j,k)/zstp/ystp)...
                    -rc(i+1,j+1,k+1)/(ystp*zstp));                   % Coefficient for vz(i,j,k)
                Lf(invy,invz-3, (muyz(i,j,k)/zstp/ystp)...
                    +rc(i,j+1,k+1)/(ystp*zstp));                    % Coefficient for vz(i-1,j,k)
            end
            
            
            % z-Stokes equation
            % Ghost vz unknowns (i=yn || j=xn) and boundary nodes (i=1, i=yn-1, j=1, j=xn-1, k=1, k=zn)
            if(i==yn || j==xn ||  k==1 || k==zn || j==1 || j==xn-1 || i==1 || i==yn-1)
                % Ghost vz unknowns (j=xn or i=yn: vz(i,j,k)=0
                if( i==yn || j==xn)
                    Lf(invz,invz, 1*kbond);                    % Coefficient for vz(i,j,k)
%                     R(invz,1)           =   0;
                end
                % Front boundary
                if((k==1) && (i<yn) && (j<xn))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet condition vz = 0
                        Lf(invz,invz, 1*kbond);
%                         R(invz,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %            Neumann codition: -vz(i,j,k+1) + vy(i,j,k) = 0;
                        Lf(invz,invz, 1*kbond);                    % Coefficient for vz(i,j,k)
                        Lf(invz,invz+xn*yn*3, -1*kbond);    %Coefficient for vz(i,j,k+1)
%                         R(invz,1)           =   0;
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
%                         R(invz,1) = 0;
                    elseif(bc==NEUMANN_BC)
                        % %           Neumann codition: vz(i,j,k) - vz(i,j,k-1) = 0;
                        Lf(invz,invz, 1*kbond);   %Coeff. for vz(i,j,k)
                        Lf(invz, invz - xn*yn*3, -1*kbond); %Coeff. for vz(i,j,k-1)
%                         R(invz,1) = 0;               %Right side
                    end
                end
                
                %             Left and right boundaries
                % Left boundary, inner points w.r.t. z-axis, (j=1, i<yn, 1<k<zn,)
                if(j==1 && i<yn && k>1 && k<zn)
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or eqv. to no slip:
                        %             % No slip vz=0: vz(i,j,k)-1/3*vz(i,j+1,k)=0
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz+yn*3, -1/3*kbond); % Coefficient for vz(i,j+1,k)
%                         R(invz,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                     elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann dvz/dx=0: vz(i,j,k)-vz(i,j+1,k)=0
                        %                 eqv. for psi boundary condtion
                        Lf(invz,invz, 1*kbond);                   % Coefficient for vz(i,j,k)
                        Lf(invz,invz+yn*3, -1*kbond);                   % Coefficient for vz(i,j+1,k)
%                         R(invz,1)=0;
                    end
                end
                % Right boundary, inner points w.r.t z-axis, (j=xn-1, i<yn, 1<k<zn)
                if(j==xn-1 && i<yn && k>1 && k<zn)
                    if(bc==DIRICHLET_BC)
                        %                 Dirichlet or no slip: vz=0: vz(i,j,k)-1/3*vz(i,j-1,k)=0
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz-yn*3, -1/3*kbond); % Coefficient for vz(i,j-1,k)
%                         R(invz,1)=0;
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc==PSI_BC)
                        % Free slip dvz/dx=0: vz(i,j,k)-vz(i,j-1,k)=0
                        %                 eqv. to psi boundary condition
                        Lf(invz,invz, 1*kbond);                    % Coefficient for vz(i,j,k)
                        Lf(invz,invz-yn*3, -1*kbond);                   % Coefficient for vz(i,j-1,k)
%                         R(invz,1)=0;
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
                        Lf(invz,invz + 3, -1/3*kbond); % Coefficient for vy(i+1,j,k)
%                         R(invz,1)          =   0;          % Right-hand-side part
                    elseif(bc==NEUMANN_BC)
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvz/dy=0: vz(i+1,j,k)-vz(i,j,k)=0
                        %                 Same for psi boundary condition!
                        Lf(invz,invz, 1*kbond);    % Coefficient for vz(i,j,k)
                        Lf(invz,invz + 3, -1*kbond);    % Coefficient for vz(i+1,j,k)
%                         R(invz,1)           =   0;          % Right-hand-side part
                    end
                end
                % Lower boundary, inner points (i=yn-1, 1<j<xn-1, 1<k<zn)
                if((i==(yn-1)) && k>1 && k<zn && j>1 && (j<(xn-1)))
                    if(bc==DIRICHLET_BC)
                        %                     Dirichlet or No slip:
                        % % No slip vz=0: vz(i,j,k)-1/3*vz(i-1,j,k)=0
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz - 3, -1/3*kbond); % Coefficient for vz(i-1,j,k)
%                         R(invz,1)         =   0; % Right part
                    elseif(bc==NEUMANN_BC )
                        %                 elseif(bc==NEUMANN_BC || bc == PSI_BC)
                        % Neumann condition dvz/dy=0: vz(i,j,k)-vz(i-1,j,k)=0
                        %                 of psi condition
                        Lf(invz,invz, 1*kbond); % Coefficient for vz(i,j,k)
                        Lf(invz,invz - 3, -1*kbond); % Coefficient for vz(i-1,j,k)
%                         R(invz,1)           =   0; % Right-hand-side part
                    end
                end
                
                %Internal nodes:
            else
                Lf(invz,invz+xn*yn*3, (2*muc(i+1,j+1,k+1)/zstp^2)...
                    +rc(i+1,j+1,k+1)/(zstp*zstp));                 % Coefficient for vz(i,j,k+1)
                Lf(invz,invz-xn*yn*3, (2*muc(i+1,j+1,k)/zstp^2)...
                    +rc(i+1,j+1,k)/(zstp*zstp));                   % Coefficient for vz(i,j,k-1)
                Lf(invz,invz, (-2*muc(i+1,j+1,k+1)/zstp^2-2*muc(i+1,j+1,k)/zstp^2)...
                    +(-muzx(i,j+1,k)/xstp^2-muzx(i,j,k)/xstp^2)...
                    +(-muyz(i+1,j,k)/ystp^2-muyz(i,j,k)/ystp^2)...
                    -rc(i+1,j+1,k+1)/(zstp*zstp) - rc(i+1,j+1,k)/(zstp*zstp)); % Coefficient for vz(i,j,k)
                
                Lf(invz,invz+yn*3, muzx(i,j+1,k)/xstp^2);                     % Coefficient for vz(i,j+1,k)
                Lf(invz,invz-yn*3, muzx(i,j,k)/xstp^2);                       % Coefficient for vz(i,j-1,k)

                Lf(invz,invz+3, muyz(i+1,j,k)/ystp^2);                     % Coefficient for vz(i+1,j,k)
                Lf(invz,invz-3, muyz(i,j,k)/ystp^2);                       % Coefficient for vz(i-1,j,k)
                

                Lf(invz,invx+yn*3, (muzx(i,j+1,k)/xstp/zstp)...
                    +rc(i+1,j+1,k+1)/(zstp*xstp));                  % Coefficient for vx(i,j+1,k)
                Lf(invz,invx+yn*3-xn*yn*3, (-muzx(i,j+1,k)/xstp/zstp)...
                    -rc(i+1,j+1,k)/(zstp*xstp));                 % Coefficient for vx(i,j+1,k-1)
                Lf(invz,invx, (-muzx(i,j,k)/xstp/zstp)...
                    -rc(i+1,j+1,k+1)/(zstp*xstp));                   % Coefficient for vx(i,j,k)
                Lf(invz,invx-xn*yn*3, (muzx(i,j,k)/xstp/zstp)...
                    +rc(i+1,j+1,k)/(zstp*xstp));                    % Coefficient for vx(i,j,k-1)
                
                Lf(invz,invy+3, (muyz(i+1,j,k)/zstp/ystp)...
                    +rc(i+1,j+1,k+1)/(zstp*ystp));                  % Coefficient for vy(i+1,j,k)
                Lf(invz,invy+3-xn*yn*3, (-muyz(i+1,j,k)/zstp/ystp)...
                    -rc(i+1,j+1,k)/(zstp*ystp));                 % Coefficient for vy(i+1,j,k-1)
                Lf(invz,invy, (-muyz(i,j,k)/zstp/ystp)...
                    -rc(i+1,j+1,k+1)/(zstp*ystp));                   % Coefficient for vy(i,j,k)
                Lf(invz,invy-xn*yn*3, (muyz(i,j,k)/zstp/ystp)...
                    +rc(i+1,j+1,k)/(zstp*ystp));                    % Coefficient for vy(i,j,k-1)
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
                invx = ((k-1)*xn*yn + (j-1)*yn + i)*3 - 2;
                invy = invx + 1;
                invz = invx + 2;
                Lf(1+NS,invx, 1*kbond); % Coefficient for vx(i,j,k)
                Lf(2+NS,invy, 1*kbond); % Coefficient for vy(i,j,k)
                Lf(3+NS,invz, 1*kbond); % Coefficient for vz(i,j,k)
            end
        end
    end
    R(1+xn*yn*3,1) = 0;
    R(2+xn*yn*3,1) = 0;
    R(3+xn*yn*3,1) = 0;
    
    
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
                invx = ((k-1)*xn*yn + (j-1)*yn+i)*3 - 2;
                invy = invx + 1;
                invz = invx + 2;
                %             For Oxy
                Lf(4+NS,invx, kbond/ystp); % Coefficient for vx(i,j,k)
                Lf(4+NS,invx-3, -kbond/ystp); % Coefficient for vx(i-1,j,k)
                Lf(4+NS,invy, -kbond/xstp); % Coefficient for vy(i,j,k)
                Lf(4+NS,invy-yn*3, kbond/xstp); % Coefficient for vy(i,j-1,k)
                
                %             For Oxz
                Lf(5+NS,invx, kbond/zstp); % Coefficient for vx(i,j,k)
                Lf(5+NS,invx-xn*yn*3, -kbond/zstp); % Coefficient for vx(i,j,k-1)
                Lf(5+NS,invz, -kbond/xstp); % Coefficient for vz(i,j,k)
                Lf(5+NS,invz-yn*3, kbond/xstp); % Coefficient for vz(i,j-1,k)
                
                %             For Oxy
                Lf(6+NS,invy, kbond/zstp); % Coefficient for vy(i,j,k)
                Lf(6+NS,invy-xn*yn*3, -kbond/zstp); % Coefficient for vy(i,j,k-1)
                Lf(6+NS,invz, -kbond/ystp); % Coefficient for vz(i,j,k)
                Lf(6+NS,invz-3, kbond/ystp); % Coefficient for vz(i-1,j,k)
                
            end
        end
    end
%     R(4+NS,1) = 0;
%     R(5+NS,1) = 0;
%     R(6+NS,1) = 0;
    
%     ecnstrnt_indx = 7;
    
end

nnz = size(i_indx,2);
display(['sparse matrix size: ' num2str(NS) '  and num of nonzeros: ' num2str(nnz)]);

% Create sparse matrix from the indices and values
L = sparse(i_indx(1:vals_cntr),j_indx(1:vals_cntr),L_vals(1:vals_cntr), ecnstrnt+NS, NS);

toc;
display('matrix built');
tic;


p = zeros(yn,xn,zn);
v = zeros(yn,xn,zn,3);
max_it = 100;
err = zeros(max_it,1);
fl = zeros(max_it,1);
for indx = 1:max_it
    rhs = setRhs(p);
%     S = L\rhs;
    [S fl(indx)] = bicgstab(L,rhs,1e-5,100);
%     [S fl(indx)] = gmres(L,rhs);
%     Rearrange the solution in matrix form
    for i=1:1:yn
        for j=1:1:xn
            for k=1:1:zn
                % Global index for vx, vy in S()
                invx = ((k-1)*xn*yn + (j-1)*yn+i)*3 - 2;
                invy = invx + 1;
                invz = invx + 2;
                % vx
                v(i,j,k,1)=S(invx);
                % vy
                v(i,j,k,2)=S(invy);
                % vz
                v(i,j,k,3)=S(invz);
            end
        end
    end
    div_v = div3D(v,xstp,ystp,zstp);
    curr_err = div_v(1:end-1,1:end-1,1:end-1)+a(2:end,2:end,2:end);
    p(2:end,2:end,2:end) = p(2:end,2:end,2:end) - rc(2:end,2:end,2:end).*(div_v(1:end-1,1:end-1,1:end-1) + a(2:end,2:end,2:end));
    err(indx) = sum(sum(sum(curr_err.^2)));
end
stem(err);
%Obtaining vector of solutions S()
% S=L\R;

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
% vy=zeros(yn,xn,zn);
% vx=zeros(yn,xn,zn);
% vz=zeros(yn,xn,zn);
% % Process all Grid points
% for i=1:1:yn
%     for j=1:1:xn
%         for k = 1:1:zn
%             % Global index for P, vx, vy and vz in S()
%             invx=((k-1)*xn*yn + (j-1)*yn+i)*3 - 2; % P
%             invy=invx+1;
%             invz=invx+2;
%             % vx
%             vx(i,j,k)=S(invx);
%             % vy
%             vy(i,j,k)=S(invy);
%             %         vz
%             vz(i,j,k)=S(invz);
%         end
%     end
% end
vx = v(:,:,:,1);
vy = v(:,:,:,2);
vz = v(:,:,:,3);

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

write_path = '';

% write_path = '/home/bkhanal/staggered3D/augLagr/real_';
% vec3DToVtk(vx,vy,vz,ax,ay,az,[write_path 'varmu3D.vtk']);
vec3DToVtk(vx1,vy1,vz1,ax,ay,az,[write_path 'varmu3D1.vtk']);
save([write_path 'results.mat'],'vx1','vy1','vz1','p','err');

end
