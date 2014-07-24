% Solution of 2D Stokes and continuity equations with finite differences
% on a regular grid using a pressure-velocity formulation
% for a medium with constant viscosity

% Clean all variables
clear;
% Clear all figures
clf;

% Numerical model parameters
% Model size, m
% xsize   =   1000000;        % Horizontal
% ysize   =   1500000;        % Vertical


% Numbers of nodes
xnum    =   15;%31;             % Horizontal
ynum    =   15;%21;             % Vertical

%bish:start
xsize   =   xnum-1;        % Horizontal
ysize   =   ynum-1;        % Vertical
%bish:end

% Grid step
xstp    =   xsize/(xnum-1); % Horizontal
ystp    =   ysize/(ynum-1); % Vertical

% Model viscosity
% eta     =   1e+21;
eta     =   1;%1e+2;

% Pressure condition in one cell (i==2 && j==3)
p0cell  =   0;

% Gravity acceleration directed downward
gy      =   0;%10; % m/s^2

% Create vectors for nodal points positions (basic nodes)
x       =   0:xstp:xsize;   % Horizontal
y       =   0:ystp:ysize;   % Vertical

% Create vectors for cell centers positions (staggered nodes)
xc      =   xstp/2:xstp:xsize-xstp/2; % Horizontal
yc      =   ystp/2:ystp:ysize-ystp/2; % Vertical

% Create array for density structure (two vertical layers)
rho     =   zeros(ynum,xnum);
for i=1:1:ynum
    for j=1:1:xnum
        % Horizontal position of the nodal point
        if(x(j)<xsize/2)
            rho(i,j)=1; %3200;  % left layer
        else
            rho(i,j)=1; %3200;  %3300;  % right layer
        end
    end
end


% Matrix of coefficients initialization
L       =   sparse(xnum*ynum*3,xnum*ynum*3);
% Vector of right part initialization
R       =   zeros(xnum*ynum*3,1);

% Computing Kcont and Kbond coefficients
kcont   =   1; %2*eta/(xstp+ystp);
kbond   =   1; %4*eta/(xstp+ystp)^2;

% Solving x-Stokes, y-Stokes and continuity equations
% x-Stokes: ETA(d2vx/dx2+d2vx/dy2)-dP/dx=0
% y-Stokes: ETA(d2vy/dx2+d2vy/dy2)-dP/dy=gy*RHO
% continuity: dvx/dx+dvy/dy=0
% Compose matrix of coefficients L()
% and vector (column) of right parts R()
% Boundary conditions: free slip
% Process all Grid points
for i=1:1:ynum
    for j=1:1:xnum
        
        % Global index for P, vx, vy in the current node
        inp     =   ((j-1)*ynum+i)*3-2; % P
        invx    =   inp+1;
        invy    =   inp+2;
        
        
        % Continuity equation
        % Ghost pressure unknowns (i=1, j=1) and boundary nodes (4 corners + one cell)
        if ( (i==1) || (j==1) || (i==2 && j==2) || (i==2 && j==xnum) || (i==ynum && j==2) || (i==ynum && j==xnum) || (i==2 && j==3))
            % Ghost pressure unknowns (i=1, j=1): P(i,j)=0
            if(i==1 || j==1)
                L(inp,inp)          =   1*kbond;    % Coefficient for P(i,j)
                R(inp,1  )          =   0;          % Right-hand-side part
            end
            % Upper and lower left corners dP/dx=0 => P(i,j)-P(i,j+1)=0
            if((i==2 && j==2) || (i==ynum && j==2))
                L(inp,inp       ) 	=   1*kbond;   % Coefficient for P(i,j)
                L(inp,inp+ynum*3)   =   -1*kbond;   % Coefficient for P(i,j+1)
                R(inp,1)            =   0;          % Right-hand-side part
            end
            % Upper and lower right corners dP/dx=0 => P(i,j)-P(i,j-1)=0
            if((i==2 && j==xnum) || (i==ynum && j==xnum))
                L(inp,inp       )	=   1*kbond;   % Coefficient for P(i,j)
                L(inp,inp-ynum*3)   =  -1*kbond;   % Coefficient for P(i,j-1)
                R(inp,1)            =   0;          % Right-hand-side part
            end
            % One cell
            if (i==2 && j==3)
                L(inp,inp)          =   1*kbond;    % Coefficient for P(i,j)
                R(inp,1)            =   p0cell;     % Right-hand-side part
            end
            %Internal nodes: dvx/dx+dvy/dy=0
        else
            %dvx/dx=(vx(i-1,j)-vx(i-1,j-1))/dx
            L(inp,invx-3       )    =    kcont/xstp; % Coefficient for vx(i-1,j)
            L(inp,invx-3-ynum*3)    =   -kcont/xstp; % Coefficient for vx(i-1,j-1)
            
            %dvy/dy=(vy(i,j-1)-vy(i-1,j-1))/dy
            L(inp,invy-ynum*3  )    =    kcont/ystp; % Coefficient for vy(i,j-1)
            L(inp,invy-3-ynum*3)    =   -kcont/ystp; % Coefficient for vy(i-1,j-1)
            % Right-hand-side part:0
            R(inp,1)                =   0;
        end
        
        % x-Stokes equation
        % Ghost vx unknowns (i=ynum) and boundary nodes (i=1, i=ynum-1, j=1, j=xnum)
        if(i==1 || i==ynum-1 || i==ynum || j==1 || j==xnum)
            % Ghost vx unknowns (i=ynum: vx(i,j)=0
            if(i==ynum)
                L(invx,invx)        =   1*kbond;    % Coefficient for vx(i,j)
                R(invx,1   )        =   0;          % Right-hand-side part
            end
            % Left and Right boundaries (j=1, j=xnum)
            if((j==1 || j==xnum) && i<ynum)
                % Free slip, No slip: vx(i,j)=0
                L(invx,invx)        =   1*kbond;    % Coefficient for vx(i,j)
                R(invx,1)           =   0;          % Right-hand-side part
            end
            % Upper boundary, inner points (i=1, 1<j<xnum)
            if(i==1 && j>1 && j<xnum)
                % Free slip dvx/dy=0: vx(i,j)-vx(i+1,j)=0
%                 L(invx,invx)        =   1*kbond;    % Coefficient for vx(i,j)
%                 L(invx,invx+3)      =  -1*kbond;    % Coefficient for vx(i+1,j)
%                 R(invx,1)           =   0;          % Right-hand-side part
%                 
                % % No slip vx=0: vx(i,j)-1/3*vx(i+1,j)=0
                L(invx,invx)       =   1*kbond;    % Coefficient for vx(i,j)
                L(invx,invx+3)     =   -1/3*kbond; % Coefficient for vx(i+1,j)
                R(invx,1)          =   10;%0;          % Right-hand-side part
            end
            % Lower boundary, iner points (i=ynum-1, 1<j<xnum)
            if(i==ynum-1 && j>1 && j<xnum)
                % Free slip dvx/dy=0: vx(i,j)-vx(i-1,j)=0
%                 L(invx,invx)        =    1*kbond; % Coefficient for vx(i,j)
%                 L(invx,invx-3)      =   -1*kbond; % Coefficient for vx(i-1,j)
%                 R(invx,1)           =   0; % Right-hand-side part
%                 % % No slip vx=0: vx(i,j)-1/3*vx(i+1,j)=0
                L(invx,invx)      =   1*kbond; % Coefficient for vx(i,j)
                L(invx,invx-3)    =   -1/3*kbond; % Coefficient for vx(i-1,j)
                R(invx,1)         =   0; % Right part
            end
            %Internal nodes: ETA(d2vx/dx2+d2vx/dy2)-dP/dx=0
        else
            %ETA(d2vx/dx2+2vx/dy2)=ETA*((vx(i,j-1)-2vx(i,j)+vx(i,j+1))/dx^2+(vx(i-1,j)-2vx(i,j)+vx(i+1,j))/dy^2)
            L(invx,invx-ynum*3)     =   eta/xstp^2;                 % Coefficient for vx(i,j-1)
            L(invx,invx-3)          =   eta/ystp^2;                 % Coefficient for vx(i-1,j)
            L(invx,invx)            =   -2*eta/xstp^2-2*eta/ystp^2; % Coefficient for vx(i,j)
            L(invx,invx+3)          =   eta/ystp^2;                 % Coefficient for vx(i+1,j)
            L(invx,invx+ynum*3)     =   eta/xstp^2;                 % Coefficient for vx(i,j+1)
            
            % -dP/dx=(P(i+1,j)-P(i+1,j+1))/dx
            L(invx,inp+3)           =   kcont/xstp;                 % Coefficient for P(i+1,j)
            L(invx,inp+3+ynum*3)    =   -kcont/xstp;                % Coefficient for P(i+1,j+1)
            % Right-hand-side part:0
            R(invx,1)               =   0;
        end
        
        % y-Stokes equation
        % Ghost vy unknowns (j=xnum) and boundary nodes (i=1, i=ynum, j=1, j=xnum-1)
        if(i==1 || i==ynum || j==1 || j==xnum-1 || j==xnum)
            % Ghost vy unknowns (j=xnum: vy(i,j)=0
            if(j==xnum)
                L(invy,invy)        =   1*kbond;                    % Coefficient for vy(i,j)
                R(invy,1)           =   0;
            end
            % Upper and lower boundaries (i=1, i=ynum)
            if((i==1 || i==ynum) && j<xnum)
                % Free slip, No slip: vy(i,j)=0
                L(invy,invy)        =   1*kbond;                    % Coefficient for vy(i,j)
                R(invy,1)           =   0;
            end
            % Left boundary, iner points (j=1, 1<i<ynum)
            if(j==1 && i>1 && i<ynum)
                % Free slip dvy/dx=0: vy(i,j)-vy(i,j+1)=0
%                 L(invy,invy)        =    1*kbond;                   % Coefficient for vy(i,j)
%                 L(invy,invy+ynum*3) =   -1*kbond;                   % Coefficient for vy(i,j+1)
%                 R(invy,1)=0;
%                 %             % No slip vy=0: vy(i,j)-1/3*vy(i,j+1)=0
                            L(invy,invy)=1*kbond; % Coefficient for vy(i,j)
                            L(invy,invy+ynum*3)=-1/3*kbond; % Coefficient for vy(i,j+1)
                            R(invy,1)=0;
            end
            % Right boundary, iner points (j=xnum-1, 1<i<ynum)
            if(j==xnum-1 && i>1 && i<ynum)
                % Free slip dvy/dx=0: vy(i,j)-vy(i,j-1)=0
%                 L(invy,invy)        =   1*kbond;                    % Coefficient for vy(i,j)
%                 L(invy,invy-ynum*3) =   -1*kbond;                   % Coefficient for vy(i,j-1)
%                 R(invy,1)=0;
                            % No slip vy=0: vy(i,j)-1/3*vy(i,j-1)=0
                            L(invy,invy)=1*kbond; % Coefficient for vy(i,j)
                            L(invy,invy-ynum*3)=-1/3*kbond; % Coefficient for vy(i,j-1)
                            R(invy,1)=0;
            end
            %Internal nodes: ETA(d2vy/dx2+d2vy/dy2)-dP/dy=-gy*RHO
        else
            %ETA(d2vx/dx2+2vx/dy2)=ETA*((vx(i,j-1)-2vx(i,j)+vx(i,j+1))/dx^2+(vx(i-1,j)-2vx(i,j)+vx(i+1,j))/dy^2)
            L(invy,invy-ynum*3)     =   eta/xstp^2;                 % Coefficient for vy(i,j-1)
            L(invy,invy-3)          =   eta/ystp^2;                 % Coefficient for vy(i-1,j)
            L(invy,invy)            =   -2*eta/xstp^2-2*eta/ystp^2; % Coefficient for vy(i,j)
            L(invy,invy+3)          =   eta/ystp^2;                 % Coefficient for vy(i+1,j)
            L(invy,invy+ynum*3)     =   eta/xstp^2;                 % Coefficient for vy(i,j+1)
            
            % -dP/dy=(P(i,j+1)-P(i+1,j+1))/dx
            L(invy,inp+ynum*3)      =   kcont/ystp;                 % Coefficient for P(i,j+1)
            L(invy,inp+3+ynum*3)    =   -kcont/ystp;                % Coefficient for P(i+1,j+1)
            
            % Right part: -RHO*gy
            R(invy,1)               =   -gy*(rho(i,j)+rho(i,j+1))/2;
        end
        
    end
end

%Obtaining vector of solutions S()
S=L\R;

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


% Compute vx,vy for internal nodes
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
figure(1);

% Plotting pressure as colormap
subplot(1,2,1);
pcolor(xc/1000,yc/1000,p(2:1:ynum,2:1:xnum)*1e-9);      % making a colormap
shading interp;     % making smooth transitions between colors
colorbar;           % showing a colorbar for the map
hold on;            % continuing plotting on the colormap
% Plotting velocity vector as arrows using internal nodes only
quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
hold off;           % stop plotting on the colormap
box on;             % making a box around the plot
title('Pressure (color,GPa), velocity (arrows)'); % title for the plot
xlabel('x, km');        % title for the horizontal axis
ylabel('y, km');        % title for the vertical axis
axis ij image ;     % directing vertical axis downward, making proper dimensions
axis([0 xsize/1000 0 ysize/1000]); % Making axes limits

% Plotting density as colormap
subplot(1,2,2);
pcolor(x/1000,y/1000,rho);      % making a colormap
shading interp;     % making smooth transitions between colors
colorbar;           % showing a colorbar for the map
hold on;            % continuing plotting on the colormap
% Plotting velocity vector as arrows using internal nodes only
quiver(x(2:1:xnum-1)/1000,y(2:1:ynum-1)/1000,vx1(2:1:ynum-1,2:1:xnum-1),vy1(2:1:ynum-1,2:1:xnum-1),'k'); % making field of arrows
hold off;           % stop plotting on the colormap
box on;             % making a box around the plot
title('Density (color, kg/m^3), velocity (arrows)');   % title for the plot
xlabel('x, km');        % title for the horizontal axis
ylabel('y, km');        % title for the vertical axis
axis ij image ;     % directing vertical axis downward, making proper dimensions
axis([0 xsize/1000 0 ysize/1000]); % Making axes limits

