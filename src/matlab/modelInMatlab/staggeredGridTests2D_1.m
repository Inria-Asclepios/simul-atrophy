% Test of staggered gridding and null spaces in 2D.
% We eliminate the Dirichlet boundary points from the matrix, by taking it
% to the right hand side.
% Regular staggered grid using a pressure-velocity formulation
% Constant viscosity case:

clear; clc;
addpath(genpath('/home/bkhanal/Documents/softwares/matlabTools/'));
addpath(genpath('/home/bkhanal/works/AdLemModel/matlab/modelInMatlab/'));

% Numbers of nodes
xn    =   23;             % Horizontal
yn    =   23;             % Vertical
N = 3*xn*yn;
% Coordinate values:
[x y] = meshgrid(1:xn,1:yn);
y = flipud(y);
% Boundary conditions:
nWallVx = 1;
nWallVy = 0;
sWallVx = 0;
sWallVy = 0;
wWallVx = 0;
wWallVy = 0;
eWallVx = 0;
eWallVy = 0;
% Call a function:
brain_mask = ones(yn,xn);
% center_mask = [round(xn/2) round(yn/2)];
% brain_mask = circleMask(center_mask,xnum/4,xnum-1,ynum-1);

% atrophy at the center of the cells
a = zeros(yn,yn);
% a = computeAtrophy(yn-1,xn-1);

% Grid step
xstp    =   1;  %xsize/(xn-1); % Horizontal
ystp    =   1;  %ysize/(yn-1); % Vertical

% Computing Kcont and Kbond coefficients
kcont   =   1;
kbond   =   1;

% Top left:
p0Cellx = 4; %indx j
p0Celly = 4; %indx i

% Process all Grid points
M = zeros(N);
R = zeros(N,1);
% maximum num of non zeros per row:
nnz = 7;
rowSt = struct('c',0,'i',0,'j',0,'k',0,'xn',xn,'yn',yn,'zn',1,'dof',3);
colSts(1) = rowSt;
for i=1:nnz
    colSts(i) = rowSt;
end
v = zeros(1,nnz);

for j=1:1:yn
    for i=1:1:xn
        rowSt.i = i;     rowSt.j = j;
        
        % x-Stokes equation: d2vx/dx^2
        rowSt.c = 0;    num = 1;
        row = getPos2D(rowSt);
        % (vx(i-1,j) + vx(i+1,j))/dx^2
        if(i>1)
            colSts(num).i = i-1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp); colSts(num).c = 0;
            num = num+1;
        end
        if(i<xn)
            colSts(num).i = i+1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp); colSts(num).c = 0;
            num = num+1;
        end
        % (vx(i,j-1) + vx(i,j+1))/dy^2
        if(j>1)
            colSts(num).i = i;  colSts(num).j = j-1;
            v(num) = 1/(ystp*ystp); colSts(num).c = 0;
            num = num+1;
        end
        if(j<yn)
            colSts(num).i = i;  colSts(num).j = j+1;
            v(num) = 1/(ystp*ystp); colSts(num).c = 0;
            num = num+1;
        end
        % -vx(i,j)/(dx^2+dy^2)
        colSts(num).i = i;  colSts(num).j = j;
        v(num) = -2/(ystp*ystp) - 2/(xstp*xstp); colSts(num).c = 0;
        num = num+1;
        
        % -dp/dx: (-p(i+1,j)+p(i,j))/xstp
        if(i<xn)
            colSts(num).i = i+1;  colSts(num).j = j;
            v(num) = -1/xstp;        colSts(num).c = 2;
            num = num+1;
        end
        colSts(num).i = i;  colSts(num).j = j;
        v(num) = 1/xstp;    colSts(num).c = 2;
        
        M(sub2ind([N N],row*ones(1,num),getPos2D(colSts(1:num)))) = v(1:num);
        
        %Update the rows which has stencils for the boundary points:
        %South wall:
        if(j==1)
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = 0;
            col = getPos2D(colSts(1));
            M(row,col) = M(row,col) - (1/ystp);
        end
        % West wall: No change in the matrix, only RHS will be updated.
        %North wall:
        if(j==yn)
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = 0;
            col = getPos2D(colSts(1));
            M(row,col) = M(row,col) + 1/(3*ystp);
        end
        % East wall: vx no change.
        %Eliminate pressure point: MUST CHECK: does this really
        % put neumann condition on the pressure ?? we simply removed 'p'
        % here!!
        if(i==xn)
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = 2;
            col = getPos2D(colSts(1));
            M(row,col) = 0;
        end
        
        % RHS: first normal interior values:
        if(i==xn)
            R(row,1) = -a(i,j)/xstp;
        else
            R(row,1) = (a(i+1,j) - a(i,j))/xstp;
        end
        %Dirichlet boundary contribution
        if(j==1) %South wall
            R(row,1) = R(row,1) - sWallVx/ystp;
        end
        if(i==1) %West wall
            R(row,1) = R(row,1) - wWallVx/xstp;
        end
        if(j==yn) %North wall
            R(row,1) = R(row,1) - nWallVx*2/(3*ystp);
        end
        if(i==xn) %East wall
            R(row,1) = R(row,1) - eWallVx/xstp;
        end
        
        % y-Stokes equation
        rowSt.c = 1;    num = 1;
        row = getPos2D(rowSt);
        % (vy(i-1,j) + vy(i+1,j))/dx^2
        if(i>1)
            colSts(num).i = i-1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp); colSts(num).c = 1;
            num = num+1;
        end
        if(i<xn)
            colSts(num).i = i+1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp);     colSts(num).c = 1;
            num = num+1;
        end
        % (vy(i,j-1) + vy(i,j+1))/dy^2
        if(j>1)
            colSts(num).i = i;  colSts(num).j = j-1;
            v(num) = 1/(ystp*ystp);     colSts(num).c = 1;
            num = num+1;
        end
        if(j<yn)
            colSts(num).i = i;  colSts(num).j = j+1;
            v(num) = 1/(ystp*ystp);     colSts(num).c = 1;
            num = num+1;
        end
        % -vy(i,j)/(dx^2+dy^2)
        colSts(num).i = i;  colSts(num).j = j;
        v(num) = -2/(ystp*ystp) - 2/(xstp*xstp);     colSts(num).c = 1;
        num = num+1;
        
        % -dp/dy: -(p(i,j+1)+p(i,j))/ystp
        if(j<yn)
            colSts(num).i = i;  colSts(num).j = j+1;
            v(num) = -1/ystp;    colSts(num).c = 2;
            num = num+1;
        end
        colSts(num).i = i;  colSts(num).j = j;
        v(num) = 1/ystp;    colSts(num).c = 2;
        
        M(sub2ind([N N],row*ones(1,num),getPos2D(colSts(1:num)))) = v(1:num);
        
        %Update the rows which has stencils for the boundary points:
        %South wall: Only RHS changes
        % West wall:
        if(i==1)
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = 1;
            col = getPos2D(colSts(1));
            M(row,col) = M(row,col) - 1/xstp;
        end
        %North wall: P gets eliminated: Check if this really makes any
        %sense or not as we are putting Neumann condition without actually
        %putting anything in the matrix!
        if(j==yn)
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = 2;
            col = getPos2D(colSts(1));
            M(row,col) = 0;
        end
        % East wall:
        if(i==xn)
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = 1;
            col = getPos2D(colSts(1));
            M(row,col) = M(row,col) + 1/(3*xstp);
        end
        
        % RHS: first normal interior values:
        if(j==yn)
            R(row,1) = -a(i,j)/ystp;
        else
            R(row,1) = (a(i,j+1) - a(i,j))/ystp;
        end
        %Dirichlet boundary contribution
        if(j==1) %South wall
            R(row,1) = R(row,1) - sWallVy/ystp;
        end
        if(i==1) %West wall
            R(row,1) = R(row,1) - 2*wWallVy/xstp;
        end
        if(j==yn) %North wall
            R(row,1) = R(row,1) - nWallVy/ystp;
        end
        if(i==xn) %East wall
            R(row,1) = R(row,1) - 2*eWallVy/(3*xstp);
        end
        
        
        %continuity equation: dvx/dx + dvy/dy =(vx(i,j)-vx(i-1,j))/dx + ...
        rowSt.c = 2;        num = 1;
        row = getPos2D(rowSt);
        %                         -vx(i-1,j)
        if(i>1)
            colSts(num).i = i-1;  colSts(num).j = j;
            v(num) = -1/xstp;      colSts(num).c = 0;
            num = num+1;
        end
        %vx(i,j)
        colSts(num).i = i;    colSts(num).j = j;
        v(num) = 1/xstp;      colSts(num).c = 0;
        num = num+1;
        
        %dvy/dy = (vy(i,j) - vy(i,j-1))/dy
        %-vy(i,j-1)
        if(j>1)
            colSts(num).i = i;  colSts(num).j = j-1;
            v(num) = -1/ystp;    colSts(num).c = 1;
            num = num+1;
        end
        %vy(i,j)
        colSts(num).i = i;    colSts(num).j = j;
        v(num) = 1/ystp;      colSts(num).c = 1;
        
        
        M(sub2ind([N N],row*ones(1,num),getPos2D(colSts(1:num)))) = v(1:num);
        
        
        %Update the rows which has stencils for the boundary points:
        %South wall: Only RHS changes.
        % West wall: Only RHS changes.
        %North wall: No change.
        % East wall: No change.
        % Fixed point:
        if(i==p0Cellx && j==p0Celly)
            M(row,:) = 0;
            M(row,row) = 1;
        end

        % RHS: first normal interior values:
        R(row,1) = -a(i,j);
        
        %Dirichlet boundary contribution
        if(j==1) %South wall
            R(row,1) = R(row,1) + sWallVy/ystp;
        end
        if(i==1) %West wall
            R(row,1) = R(row,1) + wWallVx/xstp;
        end
        %North wall: No change.
        %East wall: No change.
        
        % Fixed point:0
        if(i==p0Cellx && j==p0Celly)
            R(row,1)  =   0;
        else
            R(row,1) = -a(i,j)*kcont;
        end
    end
end

% Test number of non-zeros in each row:
% nnz = zeros(xn*yn*3); nnz(M==0) = 1; nnz = sum(~nnz,2);
% Solution
sol = M\R;
res = M*sol - R;
vx = reshape(sol(1:3:end),xn,yn);
vy = reshape(sol(2:3:end),xn,yn);
p = rot90(reshape(sol(3:3:end),xn,yn));
imagesc(p);

% Let's interpolate the values:
vxC = (vx(1:xn-1,2:yn) + vx(2:xn,2:yn))/2;
vyC = (vy(2:xn,1:yn-1) + vy(2:xn,2:yn))/2;
vxC = padarray(vxC,[1 1],0,'pre');
vyC = padarray(vyC,[1 1],0,'pre');
vxC = rot90(vxC);   vyC = rot90(vyC);
figure, quiver(x,y,vxC,vyC);
vx = rot90(vx); vy = rot90(vy);
%% Test Null Space
nullBasis = zeros(N,1);
nullBasis(3:3:end) = 1;
r = M*nullBasis;

%% Test Schur complement:
[RM fg] = reorderInterLeavedToBlockwise3Dof(M,R);
A = RM(1:2*xn*yn,1:2*xn*yn);
B = RM(1:2*xn*yn,2*xn*yn+1:end);
D = RM(2*xn*yn+1:end,1:2*xn*yn);
K = RM(2*xn*yn+1:end,2*xn*yn+1:end);
nnz = B'-D;
S = D*(A\B) - K;

[u E v] = svd(S); E = diag(E);
Kd = diag(K);
% Solution using this:
f = fg(1:2*xn*yn);
g = fg(2*xn*yn+1:end);
pS = S\((D*(A\f))-g);
uS = A\(f-B*pS);
pS = rot90(reshape(pS,xn,yn));
imagesc(pS);
vxS = reshape(uS(1:xn*yn),xn,yn);
vyS = reshape(uS(xn*yn+1:end),xn,yn);
% Let's interpolate the values:
vxCS = (vxS(1:xn-1,2:yn) + vxS(2:xn,2:yn))/2;
vyCS = (vyS(2:xn,1:yn-1) + vyS(2:xn,2:yn))/2;
vxCS = padarray(vxCS,[1 1],0,'pre');
vyCS = padarray(vyCS,[1 1],0,'pre');
vxCS = rot90(vxCS);   vyCS = rot90(vyCS);
figure, quiver(x,y,vxCS,vyCS);
vxS = rot90(vxS); vyS = rot90(vyS);
figure,
subplot(121), imagesc(vy-vyS),title('vy Difference');
subplot(122), imagesc(vx-vxS), title('vx difference');
figure,
imagesc(p-pS), title('p difference');

%% Null spaces:
[u D v] = svd(M);

currV=v(:,end);
nullP = rot90(reshape(currV(3:3:end),xn,yn));
imagesc(nullP);
nullVx = rot90(reshape(currV(1:3:end),xn,yn));
nullVy = rot90(reshape(currV(2:3:end),xn,yn));
figure, quiver(nullVx,nullVy), axis ij image;
