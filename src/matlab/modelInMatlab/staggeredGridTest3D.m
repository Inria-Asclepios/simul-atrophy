% Test of staggered gridding and null spaces in 3D.
% This is to test how using div(u) + cp = -a has an effect when c is
% non-zero in certain parts of the domain.

clear; clc;
addpath(genpath('/home/bkhanal/Documents/softwares/matlabTools/'));
addpath(genpath('/home/bkhanal/works/AdLemModel/src/matlab/modelInMatlab/'));
% Numbers of nodes
xn    =   12;             % Horizontal
yn    =   12;             % Vertical
zn    =   12;
N = 4*xn*yn*zn;

% Coordinate values:
[x y z] = meshgrid(1:xn,1:zn,1:yn);
z = flipdim(z,1);
% Boundary conditions:
nWallVx = 1;
nWallVy = 0;
nWallVz = 0;
sWallVx = 0;
sWallVy = 0;
sWallVz = 0;
wWallVx = 0;
wWallVy = 0;
wWallVz = 0;
eWallVx = 0;
eWallVy = 0;
eWallVz = 0;
fWallVx = 0;
fWallVy = 0;
fWallVz = 0;
bWallVx = 0;
bWallVy = 0;
bWallVz = 0;
% Boundary coefficient scale factors
dbScale = 1;
nullBasis = zeros(N,1);

load bMask;
% bMask = zeros(xn,yn,zn);
% bMask = ones(yn*xn,1); %let's put zeros one those points where we will
% % release the incompressibility constraint.
% bSt = createGridStencilStructure(xn,yn,1,1);
% bSt.i = 18;      bSt.j = 4;
% bMask(getPos3D(bSt)) = 0;

pMassCoeff = 10000;
% pMassCoeff = 0;
if(pMassCoeff == 0)
    bMask(bMask==0) = 1;
end

% % atrophy at the center of the cells
a = zeros(yn*xn*zn,1);
bSt = createGridStencilStructure(xn,yn,zn,1);
% 
% for k = 1:zn
%     for i=1:xn
%         for j=1:yn
%             bSt.i = i;  bSt.j = j;  bSt.k = k;
%             %         if(bMask(i,j,k)==1)
%             if(((i-4)^2 + (j-4)^2) + (k-4)^2 < 4)
%                 a(getPos3D(bSt)) = 0.2;
%             end
%         end
%     end
% end

aSt = createGridStencilStructure(xn,yn,zn,1);
aSt(2) = aSt;

% Grid step
xstp    =   1;  %xsize/(xn-1); % w-e
ystp    =   1;  %ysize/(yn-1); % f-b
zstp    =   1;  %zsize/(zn-1); % s-n
% Computing Kcont and Kbond coefficients
kcont   =   1;
kbond   =   1;

% Pin pressure point: IMP: It cannot be any of the south or west wall or
% top right corner pressure vars since they are fictious cells here!
p0Cellx = 4; %indx j
p0Celly = 4; %indx i
p0Cellz = 4; %indx k
% Process all Grid points
M = zeros(N);
R = zeros(N,1);
% maximum num of non zeros per row:
nnz = 9;
rowSt = createGridStencilStructure(xn,yn,zn,4);
colSts(1) = rowSt;
for i=1:nnz
    colSts(i) = rowSt;
end

cVx = 0;
cVy = 1;
cVz = 2;
cP = 3;
v = zeros(1,nnz);

for k = 1:1:zn
    for j=1:1:yn
        for i=1:1:xn
            rowSt.i = i;     rowSt.j = j;   rowSt.k = k;
            bSt.i = i;       bSt.j = j;     bSt.k = k;
            
            % -------------------------------------------------------------------------
            % x-Stokes equation: d2vx/dx^2
            rowSt.c = cVx;    num = 0;
            row = getPos3D(rowSt);
            % (vx(i-1,j,k) + vx(i+1,j,k))/dx^2
            if(i>1)
                num = num+1;
                colSts(num).i = i-1;  colSts(num).j = j;    colSts(num).k = k;
                v(num) = 1/(xstp*xstp); colSts(num).c = cVx;
            end
            if(i<xn)
                num = num+1;
                colSts(num).i = i+1;  colSts(num).j = j;    colSts(num).k = k;
                v(num) = 1/(xstp*xstp); colSts(num).c = cVx;
            end
            % (vx(i,j-1,k) + vx(i,j+1,k))/dy^2
            if(j>1)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j-1;    colSts(num).k = k;
                v(num) = 1/(ystp*ystp); colSts(num).c = cVx;
            end
            if(j<yn)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j+1;    colSts(num).k = k;
                v(num) = 1/(ystp*ystp); colSts(num).c = cVx;
            end
            % (vx(i,j,k-1) + vx(i,j,k+1))/dz^2
            if(k>1)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j;    colSts(num).k = k-1;
                v(num) = 1/(zstp*zstp); colSts(num).c = cVx;
            end
            if(k<zn)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j;    colSts(num).k = k+1;
                v(num) = 1/(zstp*zstp); colSts(num).c = cVx;
            end
            % -vx(i,j,k)/(dx^2+dy^2+dz^2)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j;  colSts(num).k = k;
            v(num) = -2/(zstp*zstp) -2/(ystp*ystp) - 2/(xstp*xstp);
            colSts(num).c = cVx;
            
            % -dp/dx: (-p(i+1,j+1,k+1)+p(i,j+1,K+1))/xstp
            if(k<zn)
                if(j<yn)
                    if(i<xn)
                        num = num+1;
                        colSts(num).i = i+1;  colSts(num).j = j+1;  colSts(num).k = k+1;
                        v(num) = -1/xstp;     colSts(num).c = cP;
                    end
                    num = num+1;
                    colSts(num).i = i;  colSts(num).j = j+1;    colSts(num).k = k+1;
                    v(num) = 1/xstp;    colSts(num).c = cP;
                end
            end
            
            if(num>0)
                M(sub2ind([N N],row*ones(1,num),getPos3D(colSts(1:num)))) = v(1:num);
            end
            %RESET the rows and explicitly set boundary conditions.
            %Prioritize first the ghost cells, then those that set boundary
            %condtitions without interpolation (effects orderings of the elseif
            %statements!!)
            if(j==yn || k==zn)    %Back or North wall Ghost points: vx = 0.
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                colSts(1).c = cVx;
                M(row,getPos3D(colSts(1))) = dbScale;
            elseif(i==1 || i==xn)    %West or East wall: vx = wx or ex
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                colSts(1).c = cVx;
                M(row,getPos3D(colSts(1))) = dbScale;
            elseif(j==1 || j==yn-1)    %Front or Back wall: 3vx(i,j,k) - vx(i,j+-1,k) = 2fx or 2bx
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                v(1) = 3*dbScale;   colSts(1).c = cVx;
                v(2) = -dbScale;  colSts(2).c = cVx;
                if(j==1)
                    colSts(2).i = i;    colSts(2).j = j+1;  colSts(2).k = k;
                else
                    colSts(2).i = i;    colSts(2).j = j-1;  colSts(2).k = k;
                end
                M(sub2ind([N N],row*ones(1,2),getPos3D(colSts(1:2)))) = v(1:2);
            elseif(k==1 || k==zn-1)     %South or North Wall: 3vx(i,j,k) - vx(i,j,k+-1) = 2sx or 2nx
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                v(1) = 3*dbScale;   colSts(1).c = cVx;
                v(2) = -dbScale;  colSts(2).c = cVx;
                if(k==1)
                    colSts(2).i = i;    colSts(2).j = j;  colSts(2).k = k+1;
                else
                    colSts(2).i = i;    colSts(2).j = j;  colSts(2).k = k-1;
                end
                M(sub2ind([N N],row*ones(1,2),getPos3D(colSts(1:2)))) = v(1:2);
            end
            
            % RHS: Same priority order as for the setting up of the matrix.
            if(j==yn || k==zn)    %Back and North wall, ghost point:
                R(row,1) = 0;
            elseif(i==1)    %west wall
                R(row,1) = dbScale*wWallVx;
            elseif(i==xn)    %east wall
                R(row,1) = dbScale*eWallVx;
            elseif(j==1)   %front wall
                R(row,1) = dbScale*2*fWallVx;
            elseif(j==yn-1)   %back wall
                R(row,1) = dbScale*2*bWallVx;
            elseif(k==1)   %South wall
                R(row,1) = dbScale*2*sWallVx;
            elseif(k==zn-1)   %North wall
                R(row,1) = dbScale*2*nWallVx;
            else
                aSt(1).i = i+1;   aSt(1).j = j+1;   aSt(1).k = k+1;
                aSt(2).i = i;   aSt(2).j = j+1;     aSt(2).k = k+1;
                R(row,1) = (a(getPos3D(aSt(1))) - a(getPos3D(aSt(2))))/xstp;
            end
            
            % -------------------------------------------------------------------------
            % y-Stokes equation
            rowSt.c = cVy;    num = 0;
            row = getPos3D(rowSt);
            % (vy(i-1,j) + vy(i+1,j))/dx^2
            if(i>1)
                num = num+1;
                colSts(num).i = i-1;  colSts(num).j = j;    colSts(num).k = k;
                v(num) = 1/(xstp*xstp); colSts(num).c = cVy;
            end
            if(i<xn)
                num = num+1;
                colSts(num).i = i+1;  colSts(num).j = j;    colSts(num).k = k;
                v(num) = 1/(xstp*xstp);     colSts(num).c = cVy;
            end
            % (vy(i,j-1) + vy(i,j+1))/dy^2
            if(j>1)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j-1;    colSts(num).k = k;
                v(num) = 1/(ystp*ystp);     colSts(num).c = cVy;
            end
            if(j<yn)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j+1;    colSts(num).k = k;
                v(num) = 1/(ystp*ystp);     colSts(num).c = cVy;
            end
            if(k>1)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j;    colSts(num).k = k-1;
                v(num) = 1/(zstp*zstp); colSts(num).c = cVy;
            end
            if(k<zn)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j;    colSts(num).k = k+1;
                v(num) = 1/(zstp*zstp); colSts(num).c = cVy;
            end
            % -vy(i,j)/(dx^2+dy^2)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j;  colSts(num).k = k;
            colSts(num).c = cVy;
            v(num) = -2/(zstp*zstp) -2/(ystp*ystp) - 2/(xstp*xstp);
            
            % -dp/dy: -(p(i+1,j+1)-p(i+1,j))/ystp
            if(k<zn)
                if(i<xn)
                    if(j<yn)
                        num = num+1;
                        colSts(num).i = i+1;  colSts(num).j = j+1;  colSts(num).k = k+1;
                        v(num) = -1/ystp;    colSts(num).c = cP;
                    end
                    num = num+1;
                    colSts(num).i = i+1;  colSts(num).j = j;    colSts(num).k = k+1;
                    v(num) = 1/ystp;    colSts(num).c = cP;
                end
            end
            
            if(num>0)
                M(sub2ind([N N],row*ones(1,num),getPos3D(colSts(1:num)))) = v(1:num);
            end
            %RESET the rows and explicitly set boundary conditions.
            %Prioritize first the ghost cells, then those that set boundary
            %condtitions without interpolation (effects orderings of the elseif
            %statements!!)
            if(i==xn || k==zn)        %East or North wall, ghost points:
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                colSts(1).c = cVy;
                M(row,getPos3D(colSts(1))) = dbScale;
            elseif(j==1 || j==yn)        %Front or Back wall: vy(i,j,k) = bx or fx
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                colSts(1).c = cVy;
                M(row,getPos3D(colSts(1))) = dbScale;
            elseif(i==1 || i==xn-1)    %West or East wall: 3vy(i,j,k)-vy(i+-1,j,k) = 2wy or 2ey
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                if(i==1)
                    colSts(2).i = i+1;  colSts(2).j = j;    colSts(2).k = k;
                else
                    colSts(2).i = i-1;  colSts(2).j = j;    colSts(2).k = k;
                end
                v(1) = 3*dbScale;     colSts(1).c = cVy;
                v(2) = -dbScale;     colSts(2).c = cVy;
                M(sub2ind([N N],row*ones(1,2),getPos3D(colSts(1:2)))) = v(1:2);
            elseif(k==1 || k==zn-1)    %South or North wall: 3vy(i,j,k)-vy(i,j,k+-1) = 2sy or 2ny
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                if(k==1)
                    colSts(2).i = i;    colSts(2).j = j;    colSts(2).k = k+1;
                else
                    colSts(2).i = i;  colSts(2).j = j;    colSts(2).k = k-1;
                end
                v(1) = 3*dbScale;     colSts(1).c = cVy;
                v(2) = -dbScale;     colSts(2).c = cVy;
                M(sub2ind([N N],row*ones(1,2),getPos3D(colSts(1:2)))) = v(1:2);
            end
            % RHS: Note the same order of if-elseif as in matrix rows setting.
            if(i==xn || k==zn)   %East or North wall, ghost points
                R(row,1) = 0;
            elseif(j==1)   %Front wall
                R(row,1) = dbScale*fWallVy;
            elseif(j==yn)    %Back wall
                R(row,1) = dbScale*bWallVy;
            elseif(i==1)    %west wall
                R(row,1) = dbScale*2*wWallVy;
            elseif(i==xn-1)   %east wall
                R(row,1) = dbScale*2*eWallVy;
            elseif(k==1)    %South wall
                R(row,1) = dbScale*2*sWallVy;
            elseif(k==zn-1)   %North wall
                R(row,1) = dbScale*2*nWallVy;
            else
                aSt(1).i = i+1;   aSt(1).j = j+1;   aSt(1).k = k+1;
                aSt(2).i = i+1;   aSt(2).j = j;     aSt(2).k = k+1;
                R(row,1) = (a(getPos3D(aSt(1))) - a(getPos3D(aSt(2))))/ystp;
            end
            
            % -------------------------------------------------------------------------
            % z-Stokes equation
            rowSt.c = cVz;    num = 0;
            row = getPos3D(rowSt);
            % (vz(i-1,j,k) + vy(i+1,j,k))/dx^2
            if(i>1)
                num = num+1;
                colSts(num).i = i-1;  colSts(num).j = j;    colSts(num).k = k;
                v(num) = 1/(xstp*xstp); colSts(num).c = cVz;
            end
            if(i<xn)
                num = num+1;
                colSts(num).i = i+1;  colSts(num).j = j;    colSts(num).k = k;
                v(num) = 1/(xstp*xstp);     colSts(num).c = cVz;
            end
            % (vz(i,j-1,k) + vz(i,j+1,k))/dy^2
            if(j>1)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j-1;    colSts(num).k = k;
                v(num) = 1/(ystp*ystp);     colSts(num).c = cVz;
            end
            if(j<yn)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j+1;    colSts(num).k = k;
                v(num) = 1/(ystp*ystp);     colSts(num).c = cVz;
            end
            if(k>1)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j;    colSts(num).k = k-1;
                v(num) = 1/(zstp*zstp); colSts(num).c = cVz;
            end
            if(k<zn)
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j;    colSts(num).k = k+1;
                v(num) = 1/(zstp*zstp); colSts(num).c = cVz;
            end
            % -vz(i,j)/(dx^2+dy^2+dz^2)
            num = num+1;
            colSts(num).i = i;  colSts(num).j = j;  colSts(num).k = k;
            colSts(num).c = cVz;
            v(num) = -2/(zstp*zstp) -2/(ystp*ystp) - 2/(xstp*xstp);
            
            % -dp/dz: -(p(i+1,j+1)-p(i+1,j))/ystp
            if(i<xn)
                if(j<yn)
                    if(k<zn)
                        num = num+1;
                        colSts(num).i = i+1;  colSts(num).j = j+1;  colSts(num).k = k+1;
                        v(num) = -1/zstp;    colSts(num).c = cP;
                    end
                    num = num+1;
                    colSts(num).i = i+1;  colSts(num).j = j+1;    colSts(num).k = k;
                    v(num) = 1/zstp;    colSts(num).c = cP;
                end
            end
            if(num>0)
                M(sub2ind([N N],row*ones(1,num),getPos3D(colSts(1:num)))) = v(1:num);
            end
            %RESET the rows and explicitly set boundary conditions.
            %Prioritize first the ghost cells, then those that set boundary
            %condtitions without interpolation (effects orderings of the elseif
            %statements!!)
            if(i==xn || j==yn)        %East or Back wall, ghost points:
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                colSts(1).c = cVz;
                M(row,getPos3D(colSts(1))) = dbScale;
            elseif(k==1 || k==zn)        %South or North wall: vz(i,j,k) = sz or nz
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                colSts(1).c = cVz;
                M(row,getPos3D(colSts(1))) = dbScale;
            elseif(i==1 || i==xn-1)    %West or East wall: 3vz(i,j,k)-vz(i+-1,j,k) = 2wz or 2ez
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                if(i==1)
                    colSts(2).i = i+1;  colSts(2).j = j;    colSts(2).k = k;
                else
                    colSts(2).i = i-1;  colSts(2).j = j;    colSts(2).k = k;
                end
                v(1) = 3*dbScale;     colSts(1).c = cVz;
                v(2) = -dbScale;     colSts(2).c = cVz;
                M(sub2ind([N N],row*ones(1,2),getPos3D(colSts(1:2)))) = v(1:2);
            elseif(j==1 || j==yn-1)    %Front or Back wall: 3vz(i,j,k)-vz(i,j+-1,k) = 2fz or 2bz
                M(row,:) = 0;
                colSts(1).i = i;    colSts(1).j = j;    colSts(1).k = k;
                if(j==1)
                    colSts(2).i = i;  colSts(2).j = j+1;  colSts(2).k = k;
                else
                    colSts(2).i = i;  colSts(2).j = j-1;    colSts(2).k = k;
                end
                v(1) = 3*dbScale;     colSts(1).c = cVz;
                v(2) = -dbScale;     colSts(2).c = cVz;
                M(sub2ind([N N],row*ones(1,2),getPos3D(colSts(1:2)))) = v(1:2);
            end
            % RHS: Note the same order of if-elseif as in matrix rows setting.
            if(i==xn || j==yn)   %East or Back wall, ghost points
                R(row,1) = 0;
            elseif(k==1)   %South wall
                R(row,1) = dbScale*sWallVy;
            elseif(k==zn)    %North wall
                R(row,1) = dbScale*nWallVy;
            elseif(i==1)    %west wall
                R(row,1) = dbScale*2*wWallVy;
            elseif(i==xn-1)   %east wall
                R(row,1) = dbScale*2*eWallVy;
            elseif(j==1)    %Front wall
                R(row,1) = dbScale*2*fWallVy;
            elseif(j==yn-1)  %Back wall
                R(row,1) = dbScale*2*bWallVy;
            else
                aSt(1).i = i+1;   aSt(1).j = j+1;   aSt(1).k = k+1;
                aSt(2).i = i+1;   aSt(2).j = j+1;     aSt(2).k = k;
                R(row,1) = (a(getPos3D(aSt(1))) - a(getPos3D(aSt(2))))/zstp;
            end
            
            % ------------------------------------------------------------------------
            %continuity equation: dvx/dx + dvy/dy + dvz/dz =(vx(i,j-1)-vx(i-1,j-1))/dx + ...
            rowSt.c = cP;        num = 0;
            row = getPos3D(rowSt);
            %                         -vx(i-1,j-1,k-1)
            if(k>1)
                if(j>1)
                    if(i>1)
                        num = num+1;
                        colSts(num).i = i-1;  colSts(num).j = j-1;  colSts(num).k = k-1;
                        v(num) = -1/xstp;      colSts(num).c = cVx;
                    end
                    %vx(i,j-1,k-1)
                    num = num+1;
                    colSts(num).i = i;    colSts(num).j = j-1;  colSts(num).k = k-1;
                    v(num) = 1/xstp;      colSts(num).c = cVx;
                end
            end
            
            %dvy/dy = (vy(i-1,j,k-1) - vy(i-1,j-1,k-1))/dy
            %-vy(i-1,j-1,k-1)
            if(k>1)
                if(i>1)
                    if(j>1)
                        num = num+1;
                        colSts(num).i = i-1;  colSts(num).j = j-1;  colSts(num).k = k-1;
                        v(num) = -1/ystp;    colSts(num).c = cVy;
                    end
                    %vy(i-1,j,k-1)
                    num = num+1;
                    colSts(num).i = i-1;    colSts(num).j = j;      colSts(num).k = k-1;
                    v(num) = 1/ystp;      colSts(num).c = cVy;
                end
            end
            
            %dvz/dz = (vz(i-1,j-1,k) - vz(i-1,j-1,k-1))/dz
            %-vz(i-1,j-1,k-1)
            if(j>1)
                if(i>1)
                    if(k>1)
                        num = num+1;
                        colSts(num).i = i-1;  colSts(num).j = j-1;  colSts(num).k = k-1;
                        v(num) = -1/zstp;    colSts(num).c = cVz;
                    end
                    %vz(i-1,j-1,k)
                    num = num+1;
                    colSts(num).i = i-1;    colSts(num).j = j-1;      colSts(num).k = k;
                    v(num) = 1/zstp;      colSts(num).c = cVz;
                end
            end
            
            % Release incompressibility at selected points by adding
            % coefficient for pressure
            if(~bMask(getPos3D(bSt)))
                num = num+1;
                colSts(num).i = i;  colSts(num).j = j;  colSts(num).k = k;
                v(num) = pMassCoeff;    colSts(num).c = cP;
            end
            
            if(num>0)
                M(sub2ind([N N],row*ones(1,num),getPos3D(colSts(1:num)))) = v(1:num);
            end
            
            % Boundary conditions:
            % The ghost points are those that never come into play and in being
            % linked to other parts of the equation. In this sense The first
            % row,col are all ghost pressure points in addition to the 4
            % corners. But let's see what effect we will have in letting a
            % neumann or dirichlet condtion to the 4 corners:
            
            if(k==1 || j==1 || i==1)    %South or Front or West wall: ghost points
                M(row,:) = 0;
                M(row,row) = dbScale;
            elseif(i==2 && (j==2 || j==yn || k==2 || k==zn)) % west corners
                M(row,:) = 0;
                M(row,row) = dbScale;
            elseif(i==xn && (j==2 || j==yn || k==2 || k==zn)) %east corners
                M(row,:) = 0;
                M(row,row) = dbScale; %For Dirichlet condition.
            elseif(j==2 && (k==2 || k==zn)) %front wall horizontal corners
                M(row,:) = 0;
                M(row,row) = dbScale;
            elseif(j==yn && (k==2 || k==zn)) %back wall horizontal corners
                M(row,:) = 0;
                M(row,row) = dbScale;
%             elseif(i==p0Cellx && j==p0Celly && k==p0Cellz) %Fixed Point
%                 if(pMassCoeff==0)            %Pin only if pressure mass on diagonal not used.
%                     M(row,:) = 0;
%                     M(row,row) = dbScale;
%                 else                        %If used, still let's see how null space looks.
%                     if(bMask(getPos3D(bSt)))
%                         nullBasis(row)=1;
%                     end
%                 end
            elseif(bMask(getPos3D(bSt))) %Else here just to update the null space basis.
                nullBasis(row) = 1;
            end
            
            % RHS
            if(k==1 || i==1 || j==1)    %south or Front or West walls: ghost points
                R(row,1) = 0;
            elseif(i==2 && (j==2 || j== yn || k==2 || k==zn)) %west corners
                R(row,1) = 0;
            elseif(i==xn && (j==2 || j== yn || k==2 || k==zn)) %east corners
                R(row,1) = 0;
            elseif(j==2 && (k==2 || k==zn)) %front wall horizontal corners
                R(row,1) = 0;
            elseif(j==yn && (k==2 || k==zn)) %back wall horizontal corners
                R(row,1) = 0;
%             elseif(i==p0Cellx && j==p0Celly && k==p0Cellz)         % Fixed point:
%                 if(pMassCoeff==0)
%                     R(row,1) = 0;
%                 else
%                     aSt(1).i = i;   aSt(1).j = j;   aSt(1).k = k;
%                     R(row,1) = -a(getPos3D(aSt(1)));
%                 end
            else
                aSt(1).i = i;   aSt(1).j = j;       aSt(1).k = k;
                R(row,1) = -a(getPos3D(aSt(1)));
            end
        end
    end
end
% rank(M)
% Test number of non-zeros in each row:
% nnz = zeros(xn*yn*3); nnz(M==0) = 1; nnz = sum(~nnz,2);
%% Solution
sol = M\R;
res = M*sol - R;
vx = reshape(sol(cVx+1:4:end),xn,yn,zn);
vy = reshape(sol(cVy+1:4:end),xn,yn,zn);
vz = reshape(sol(cVz+1:4:end),xn,yn,zn);
p = reshape(sol(cP+1:4:end),xn,yn,zn);

vec3DToVtk(vx,vy,vz,x,y,z,'velFile.vtk');
% % Let's interpolate the values:
% vxC = (vx(1:xn-1,2:yn) + vx(2:xn,2:yn))/2;
% vyC = (vy(2:xn,1:yn-1) + vy(2:xn,2:yn))/2;
% vxC = padarray(vxC,[1 1],0,'pre');
% vyC = padarray(vyC,[1 1],0,'pre');
% vxC = rot90(vxC);   vyC = rot90(vyC);
% figure, quiver(x,y,vxC,vyC);
% vx = rot90(vx); vy = rot90(vy);
% quiver(x,y,vx,vy);
%% Test Null Space
r = M*nullBasis;
tst = find(r~=0);

%% Test Schur complement:
[RM fg]= reorderInterLeavedToBlockwise4Dof(M,R);
A = RM(1:3*xn*yn*zn,1:3*xn*yn*zn);
B = RM(1:3*xn*yn*zn,3*xn*yn*zn+1:end);
D = RM(3*xn*yn*zn+1:end,1:3*xn*yn*zn);
K = RM(3*xn*yn*zn+1:end,3*xn*yn*zn+1:end);
nnz = B'-D;
S = D*(A\B) - K;

[u E v] = svd(S); E = diag(E);
Kd = diag(K);
%% Let's solve using this:
f = fg(1:3*xn*yn*zn);
g = fg(3*xn*yn*zn+1:end);
pS = S\((D*(A\f))-g);
uS = A\(f-B*pS);
pS = (reshape(pS,xn,yn,zn));

vxS = reshape(uS(1:xn*yn*zn),xn,yn,zn);
vyS = reshape(uS(xn*yn*zn+1:2*xn*yn*zn),xn,yn,zn);
vzS = reshape(uS(2*xn*yn*zn+1:end),xn,yn,zn);
% % Let's interpolate the values:
% vxCS = (vxS(1:xn-1,2:yn) + vxS(2:xn,2:yn))/2;
% vyCS = (vyS(2:xn,1:yn-1) + vyS(2:xn,2:yn))/2;
% vxCS = padarray(vxCS,[1 1],0,'pre');
% vyCS = padarray(vyCS,[1 1],0,'pre');
% vxCS = rot90(vxCS);   vyCS = rot90(vyCS);
% figure, quiver(x,y,vxCS,vyCS);
% vxS = rot90(vxS); vyS = rot90(vyS);
% figure,
% subplot(121), imagesc(vy-vyS),title('vy Difference');
% subplot(122), imagesc(vx-vxS), title('vx difference');
% figure,
% imagesc(p-pS), title('p difference');

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
        pos = getPos3D(posSt);
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
        
        divVal(i,j) = sol(getPos3D(vxPos(1))) + sol(getPos3D(vyPos(1)));
        if(i>1)
            vxPos(2).i = i-1; vxPos(2).j = j;
            divVal(i,j) = divVal(i,j) - sol(getPos3D(vxPos(2)));
        end
        if(j>1)
            vyPos(2).i = i; vyPos(2).j = j-1;
            divVal(i,j) = divVal(i,j) - sol(getPos3D(vyPos(2)));
        end
    end
end
divVal = rot90(divVal);
figure, subplot(121), imagesc(divVal), subplot(122), imagesc(rot90(reshape(a,xn,yn)));