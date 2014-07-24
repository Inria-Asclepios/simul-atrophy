% Test of staggered gridding and null spaces in 2D.
% This is to test finally the Taras Method! 
% else. Only difference with Taras is the co-ordinate system I use here,
% which is just a 90 degree rotation.
% Regular staggered grid using a pressure-velocity formulation
% Constant viscosity case:

clear; clc;
addpath(genpath('/home/bkhanal/Documents/softwares/matlabTools/'));
addpath(genpath('/home/bkhanal/works/AdLemModel/matlab/modelInMatlab/'));

% Numbers of nodes
% load bMask;
load btest2d.mat;
[xn yn] = size(btest2d);
% xn    =   24;             % Horizontal
% yn    =   24;             % Vertical
N = 3*xn*yn;

% Coordinate values:
[x y] = meshgrid(1:xn,1:yn);
y = flipud(y);
% Boundary conditions:
nWallVx = 0;
nWallVy = 0;
sWallVx = 0;
sWallVy = 0;
wWallVx = 0;
wWallVy = 0;
eWallVx = 0;
eWallVy = 0;
% Boundary coefficient scale factors
dbScale = 1;
% Call a function:
nullBasis = zeros(N,1);


% bMask = ones(yn*xn,1); %let's put zeros one those points where we will 
% % release the incompressibility constraint.
% bSt = createGridStencilStructure(xn,yn,1,1);
% bSt.i = 18;      bSt.j = 4;
% bMask(getPos2D(bSt)) = 0;

pMassCoeff = 1;  
bMask(bMask==0) = 2;  %Use this if you want all non-CSF area to conserve the density.
% pMassCoeff = 0;

% pMassCoeff = 0;
if(pMassCoeff == 0)
    bMask(bMask<2) = 2; %2 => brain parenchyma, less than that are CSF and non-brain areas.
end

% a = a1(:);
% a = a2(:);
% a = a3(:);
a = a4(:);

% %% plot
% climMin = min([min(a1(:)),min(a2(:)),min(a3(:)),min(a4(:))]);
% climMax = max([max(a1(:)),max(a2(:)),max(a3(:)),max(a4(:))]);
% clim = [climMin climMax];
% subplot(221), imagesc(rot90(a1), clim), axis image;
% subplot(222), imagesc(rot90(a2), clim), axis image;
% subplot(223), imagesc(rot90(a3), clim), axis image;
% subplot(224), imagesc(rot90(a4), clim), axis image;

bSt = createGridStencilStructure(xn,yn,1,1);
aSt = createGridStencilStructure(xn,yn,1,1);
aSt(2) = aSt;

% Grid step
xstp    =   1;  %xsize/(xn-1); % Horizontal
ystp    =   1;  %ysize/(yn-1); % Vertical

% Computing Kcont and Kbond coefficients
kcont   =   1;
kbond   =   1;

% Pin pressure point: IMP: It cannot be any of the south or west wall or
% top right corner pressure vars since they are fictious cells here!
p0Cellx = 3; %indx j
p0Celly = 2; %indx i

% Process all Grid points
M = zeros(N);
R = zeros(N,1);
% maximum num of non zeros per row:
nnz = 7;
rowSt = createGridStencilStructure(xn,yn,1,3);
colSts(1) = rowSt;
for i=1:nnz
    colSts(i) = rowSt;
end

cVx = 0;
cVy = 1;
cP = 2;
v = zeros(1,nnz);

for j=1:1:yn
    for i=1:1:xn
        rowSt.i = i;     rowSt.j = j;
        bSt.i = i;       bSt.j = j;
        
        % -------------------------------------------------------------------------
        % x-Stokes equation: d2vx/dx^2
        rowSt.c = cVx;    num = 0;
        row = getPos2D(rowSt);
        % (vx(i-1,j) + vx(i+1,j))/dx^2
        if(i>1)
            num = num+1;
            colSts(num).i = i-1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp); colSts(num).c = cVx;
        end
        if(i<xn)
            num = num+1;
            colSts(num).i = i+1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp); colSts(num).c = cVx;
        end
        % (vx(i,j-1) + vx(i,j+1))/dy^2
        if(j>1)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j-1;
            v(num) = 1/(ystp*ystp); colSts(num).c = cVx;
        end
        if(j<yn)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j+1;
            v(num) = 1/(ystp*ystp); colSts(num).c = cVx;
        end
        % -vx(i,j)/(dx^2+dy^2)
        num = num+1;
        colSts(num).i = i;  colSts(num).j = j;
        v(num) = -2/(ystp*ystp) - 2/(xstp*xstp); colSts(num).c = cVx;
        
        % -dp/dx: (-p(i+1,j+1)+p(i,j+1))/xstp
        if(j<yn)
            if(i<xn)
                num = num+1;
                colSts(num).i = i+1;  colSts(num).j = j+1;
                v(num) = -1/xstp;     colSts(num).c = cP;
            end
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j+1;
            v(num) = 1/xstp;    colSts(num).c = cP;
        end
        
        if(num>0)
            M(sub2ind([N N],row*ones(1,num),getPos2D(colSts(1:num)))) = v(1:num);
        end
        %RESET the rows and explicitly set boundary conditions.
        %Prioritize first the ghost cells, then those that set boundary
        %condtitions without interpolation (effects orderings of the elseif
        %statements!!)
        if(j==yn)    %North wall Ghost points: vx = 0.
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = cVx;
            M(row,getPos2D(colSts(1))) = dbScale;
        elseif(i==1)    %West wall: vx = wx
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = cVx;
            M(row,getPos2D(colSts(1))) = dbScale;
        elseif(i==xn)   %East wall: vx = ex
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = cVx;
            M(row,getPos2D(colSts(1))) = dbScale;
        elseif(j==1)    %South wall: 3vx(i,j) - vx(i,j+1) = 2sx
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;
%             v(1) = 3*dbScale;   colSts(1).c = cVx;
            v(1) = dbScale;   colSts(1).c = cVx;
            colSts(2).i = i;    colSts(2).j = j+1;
%             v(2) = -dbScale;  colSts(2).c = cVx;
            v(2) = -dbScale/3;  colSts(2).c = cVx;
            M(sub2ind([N N],row*ones(1,2),getPos2D(colSts(1:2)))) = v(1:2);
        elseif(j==yn-1)   %North wall: 3vx(i,j) - vx(i,j-1) = 2nx
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;
%             v(1) = 3*dbScale;   colSts(1).c = cVx;
            v(1) = dbScale;   colSts(1).c = cVx;
            colSts(2).i = i;    colSts(2).j = j-1;
%             v(2) = -1*dbScale;  colSts(2).c = cVx;
            v(2) = -1*dbScale/3;  colSts(2).c = cVx;
            M(sub2ind([N N],row*ones(1,2),getPos2D(colSts(1:2)))) = v(1:2);
        end
        
        % RHS: Same priority order as for the setting up of the matrix.
        if(j==yn)    %North wall, ghost point:
            R(row,1) = 0;
        elseif(i==1)    %west wall
            R(row,1) = dbScale*wWallVx;
        elseif(i==xn)    %east wall
            R(row,1) = dbScale*eWallVx;
        elseif(j==1)   %south wall
            R(row,1) = dbScale*2*sWallVx;
        elseif(j==yn-1)   %north wall
            R(row,1) = dbScale*2*nWallVx;
        else
            aSt(1).i = i+1;   aSt(1).j = j+1;
            aSt(2).i = i;   aSt(2).j = j+1;
            R(row,1) = (a(getPos2D(aSt(1))) - a(getPos2D(aSt(2))))/xstp;
        end
        % -------------------------------------------------------------------------
        % y-Stokes equation
        rowSt.c = cVy;    num = 0;
        row = getPos2D(rowSt);
        % (vy(i-1,j) + vy(i+1,j))/dx^2
        if(i>1)
            num = num+1;
            colSts(num).i = i-1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp); colSts(num).c = cVy;
        end
        if(i<xn)
            num = num+1;
            colSts(num).i = i+1;  colSts(num).j = j;
            v(num) = 1/(xstp*xstp);     colSts(num).c = cVy;
        end
        % (vy(i,j-1) + vy(i,j+1))/dy^2
        if(j>1)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j-1;
            v(num) = 1/(ystp*ystp);     colSts(num).c = cVy;
        end
        if(j<yn)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j+1;
            v(num) = 1/(ystp*ystp);     colSts(num).c = cVy;
        end
        % -vy(i,j)/(dx^2+dy^2)
        num = num+1;
        colSts(num).i = i;  colSts(num).j = j;
        v(num) = -2/(ystp*ystp) - 2/(xstp*xstp);     colSts(num).c = cVy;
        
        % -dp/dy: -(p(i+1,j+1)-p(i+1,j))/ystp
        if(i<xn)
            if(j<yn)
                num = num+1;
                colSts(num).i = i+1;  colSts(num).j = j+1;
                v(num) = -1/ystp;    colSts(num).c = cP;
            end
            num = num+1;
            colSts(num).i = i+1;  colSts(num).j = j;
            v(num) = 1/ystp;    colSts(num).c = cP;
        end
        
        if(num>0)
            M(sub2ind([N N],row*ones(1,num),getPos2D(colSts(1:num)))) = v(1:num);
        end
        %RESET the rows and explicitly set boundary conditions.
        %Prioritize first the ghost cells, then those that set boundary
        %condtitions without interpolation (effects orderings of the elseif
        %statements!!)
        if(i==xn)        %East wall, ghost points:
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = cVy;
            M(row,getPos2D(colSts(1))) = dbScale;
        elseif(j==1)        %South wall: vy(i,j) = sx
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = cVy;
            M(row,getPos2D(colSts(1))) = dbScale;
        elseif(j==yn)   %North wall: vy(i,j) = ny
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;    colSts(1).c = cVy;
            M(row,getPos2D(colSts(1))) = dbScale;
        elseif(i==1)    %West wall: 3vy(i,j)-vy(i+1,j) = 2wy
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;
%             v(1) = 3*dbScale;     colSts(1).c = cVy;
            v(1) = dbScale;     colSts(1).c = cVy;
            colSts(2).i = i+1;  colSts(2).j = j;
%             v(2) = -dbScale;     colSts(2).c = cVy;
            v(2) = -dbScale/3;     colSts(2).c = cVy;
            M(sub2ind([N N],row*ones(1,2),getPos2D(colSts(1:2)))) = v(1:2);
        elseif(i==xn-1)   %East wall: 3vy(i,j) - vy(i-1,j) = 2ey
            M(row,:) = 0;
            colSts(1).i = i;    colSts(1).j = j;
%             v(1) = 3*dbScale;   colSts(1).c = cVy;
            v(1) = dbScale;   colSts(1).c = cVy;
            colSts(2).i = i-1;    colSts(2).j = j;
%             v(2) = -1*dbScale;  colSts(2).c = cVy;
            v(2) = -1*dbScale/3;  colSts(2).c = cVy;
            M(sub2ind([N N],row*ones(1,2),getPos2D(colSts(1:2)))) = v(1:2);
        end
        
        % RHS: Note the same order of if-elseif as in matrix rows setting.
        if(i==xn)   %East wall, ghost points
            R(row,1) = 0;
        elseif(j==1)   %south wall
            R(row,1) = dbScale*sWallVy;
        elseif(j==yn)    %north wall
            R(row,1) = dbScale*nWallVy;
        elseif(i==1)    %west wall
            R(row,1) = dbScale*2*wWallVy;
        elseif(i==xn-1)   %east wall
            R(row,1) = dbScale*2*eWallVy;
        else
            aSt(1).i = i+1;   aSt(1).j = j+1;
            aSt(2).i = i+1;   aSt(2).j = j;
            R(row,1) = (a(getPos2D(aSt(1))) - a(getPos2D(aSt(2))))/ystp;
        end
        
        % ------------------------------------------------------------------------
        %continuity equation: dvx/dx + dvy/dy =(vx(i,j-1)-vx(i-1,j-1))/dx + ...
        rowSt.c = cP;        num = 0;
        row = getPos2D(rowSt);
        %                         -vx(i-1,j-1)
        if(j>1)
            if(i>1)
                num = num+1;
                colSts(num).i = i-1;  colSts(num).j = j-1;
                v(num) = -1/xstp;      colSts(num).c = cVx;
            end
            %vx(i,j-1)
            num = num+1;
            colSts(num).i = i;    colSts(num).j = j-1;
            v(num) = 1/xstp;      colSts(num).c = cVx;
        end
        
        %dvy/dy = (vy(i-1,j) - vy(i-1,j-1))/dy
        %-vy(i-1,j-1)
        if(i>1)
            if(j>1)
                num = num+1;
                colSts(num).i = i-1;  colSts(num).j = j-1;
                v(num) = -1/ystp;    colSts(num).c = cVy;
            end
            %vy(i-1,j)
            num = num+1;
            colSts(num).i = i-1;    colSts(num).j = j;
            v(num) = 1/ystp;      colSts(num).c = cVy;
        end
        
        % Release incompressibility at selected points by adding 
        % coefficient for pressure
        if(bMask(getPos2D(bSt))<2)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j;
            v(num) = pMassCoeff;    colSts(num).c = cP;
        end
        
        if(num>0)
            M(sub2ind([N N],row*ones(1,num),getPos2D(colSts(1:num)))) = v(1:num);
        end
        
        % Boundary conditions:
        % The ghost points are those that never come into play and in being
        % linked to other parts of the equation. In this sense The first
        % row,col are all ghost pressure points in addition to the 4
        % corners. But let's see what effect we will have in letting a
        % neumann or dirichlet condtion to the 4 corners:
        if(j==1 || i==1)    %South wall or west wall: ghost points
            M(row,:) = 0;
            M(row,row) = dbScale;
        elseif(i==2 && (j==2 || j==yn)) % west corners
            M(row,:) = 0;
            M(row,row) = dbScale;  %For Dirichlet condition.
        elseif(i==xn && (j==2 || j==yn)) %east corners
            M(row,:) = 0;
            M(row,row) = dbScale; %For Dirichlet condition.
        elseif(i==p0Cellx && j==p0Celly && pMassCoeff==0) %Fixed Point
        %Pin only if pressure mass on diagonal not used.
                M(row,:) = 0;
                M(row,row) = dbScale;
        elseif(pMassCoeff~=0) %If no pressure mass coeff, then there is the constant pressure null space.
            nullBasis(row) = 1;
        end
        
        % RHS
        if(i==1 || j==1)    %south and west walls: ghost points
            R(row,1) = 0;
        elseif(i==2 && (j==2 || j== yn)) %west corners
            R(row,1) = 0;
        elseif(i==xn && (j==2 || j== yn)) %east corners
            R(row,1) = 0;
        elseif(i==p0Cellx && j==p0Celly && pMassCoeff==0) % Fixed point:
            R(row,1) = 0;
        else
            aSt(1).i = i;   aSt(1).j = j;
            R(row,1) = -a(getPos2D(aSt(1)));
        end
    end
end

% rank(M)
% Test number of non-zeros in each row:
% nnz = zeros(xn*yn*3); nnz(M==0) = 1; nnz = sum(~nnz,2);
% Solution
sol = M\R;
res = M*sol - R;
vx = reshape(sol(cVx+1:3:end),xn,yn);
vy = reshape(sol(cVy+1:3:end),xn,yn);
p = rot90(reshape(sol(cP+1:3:end),xn,yn));
% subplot(121), imagesc(p), axis image, title('pressure map');
imagesc(rot90(reshape(a,xn,yn))), axis image, title('atrophy map'), colormap(gray);

% Let's interpolate the values:
vxC = (vx(1:xn-1,2:yn) + vx(2:xn,2:yn))/2;
vyC = (vy(2:xn,1:yn-1) + vy(2:xn,2:yn))/2;
vxC = padarray(vxC,[1 1],0,'pre');
vyC = padarray(vyC,[1 1],0,'pre');
vxC = rot90(vxC);   vyC = rot90(vyC);
figure, quiver(x,y,vxC,vyC), axis image, title('velocity field');
vx = rot90(vx); vy = rot90(vy);
% quiver(x,y,vx,vy);
%% Test Null Space
r = M*nullBasis;
tst = find(r~=0);

%% Test Schur complement:
[RM fg]= reorderInterLeavedToBlockwise3Dof(M,R);
A = RM(1:2*xn*yn,1:2*xn*yn);
B = RM(1:2*xn*yn,2*xn*yn+1:end);
D = RM(2*xn*yn+1:end,1:2*xn*yn);
K = RM(2*xn*yn+1:end,2*xn*yn+1:end);
nnz = B'-D;
S = D*(A\B) - K;

[u E v] = svd(S); E = diag(E);
Kd = diag(K);
%Let's solve using this:
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

currV=v(:,end-1);
nullP = rot90(reshape(currV(cP+1:3:end),xn,yn));
imagesc(nullP);
nullVx = rot90(reshape(currV(cVx+1:3:end),xn,yn));
nullVy = rot90(reshape(currV(cVy+1:3:end),xn,yn));
figure, quiver(x,y,nullVx,nullVy);

%% Test Quiver plots
vxTest = zeros(yn,xn); vyTest = zeros(yn,xn);
posSt.xn = xn; posSt.yn = yn; posSt.dof = 1;
for j=1:yn
    for i = 1:xn
        posSt.i = i;  posSt.j = j;  posSt.c = 0;
        pos = getPos2D(posSt);
        vxTest(i,j) = i;
        vyTest(i,j) = j;
    end
end
vxTest = rot90(vxTest); vyTest = rot90(vyTest);
quiver(x,y,vxTest,vyTest);

%% Test divergence of the solution
divVal = zeros(yn,xn);
for i=1:2
    vxPos(i) = createGridStencilStructure(xn,yn,1,3);
    vyPos(i) = createGridStencilStructure(xn,yn,1,3);
end
for i = 1:2
    vxPos(i).c = cVx;
    vyPos(i).c = cVy;
end
for j=1:yn
    for i = 1:xn
        vxPos(1).i = i;  vxPos(1).j = j;
        vyPos(1).i = i;  vyPos(1).j = j;
        
        divVal(i,j) = sol(getPos2D(vxPos(1))) + sol(getPos2D(vyPos(1)));
        if(i>1)
            vxPos(2).i = i-1; vxPos(2).j = j;
            divVal(i,j) = divVal(i,j) - sol(getPos2D(vxPos(2)));
        end
        if(j>1)
            vyPos(2).i = i; vyPos(2).j = j-1;
            divVal(i,j) = divVal(i,j) - sol(getPos2D(vyPos(2)));
        end
    end
end
divVal = rot90(divVal);
figure, subplot(121), imagesc(divVal), subplot(122), imagesc(rot90(reshape(a,xn,yn)));