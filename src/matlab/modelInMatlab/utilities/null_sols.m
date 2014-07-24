LL = full(L);
[ax bx cx] = svd(LL);
% null_sol = cx(:,end-1);
% null_sol = cx(:,1249);
null_sol = cx(:,end);
p=zeros(yn,xn);
vy=zeros(yn,xn);
vx=zeros(yn,xn);
% Process all Grid points
for i=1:1:yn
    for j=1:1:xn
        % Global index p, vx, vy in S()
%         inp=((j-1)*yn+i)*3-2; % P
%         invx=inp+1;
%         invy=inp+2;
% P
%         p(i,j)=null_sol(inp)*kcont;
        % vx
%         vx(i,j)=null_sol(invx);
        % vy
%         vy(i,j)=null_sol(invy);
        
        iinvx=((j-1)*yn+i)*2-1; % vx
        iinvy = iinvx+1;
        vx(i,j) = null_sol(invx);
        vy(i,j)=null_sol(invy);
    end
end
% Compute vx,vy for internal nodes
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


%Plotting solution
% Making new figure
figure(1);

% % Plotting pressure as colormap
% % subplot(1,2,1);
% pcolor(xc/1000,yc/1000,p(2:1:yn,2:1:xn)*1e-9);      % making a colormap
% shading interp;     % making smooth transitions between colors
% colorbar;           % showing a colorbar for the map
% hold on;            % continuing plotting on the colormap
% Plotting velocity vector as arrows using internal nodes only
quiver(x(2:1:xn-1)/1000,y(2:1:yn-1)/1000,vx1(2:1:yn-1,2:1:xn-1),vy1(2:1:yn-1,2:1:xn-1),'k'); % making field of arrows
hold off;           % stop plotting on the colormap
% box on;             % making a box around the plot
% title('Pressure (color,GPa), velocity (arrows)'); % title for the plot
% xlabel('x, m');        % title for the horizontal axis
% ylabel('y, m');        % title for the vertical axis
axis ij image ;     % directing vertical axis downward, making proper dimensions
% axis([0 xsize/1000 0 ysize/1000]); % Making axes limits

% soltn(:,:,1) = vx;
% soltn(:,:,2) = vy;
% div_ans = div2D_new(soltn,1,1);
% figure, imagesc(div_ans), title('divergence of displ. field corresponding to null space');

% Plotting singular values and the ratio of the adjacent singular values:
% sing_vals = diag(bx);
% figure,
% subplot(121), stem(sing_vals), title('singular values');
% sing_vals(sing_vals==0) = inf; %to enforce those divided by zero to be zero.
% sing_vals_ratio = sing_vals(1:end-1)/sing_vals(2:end);
% subplot(122), stem(sing_vals_ratio), title('ratios of adjacent singular values');
% 

