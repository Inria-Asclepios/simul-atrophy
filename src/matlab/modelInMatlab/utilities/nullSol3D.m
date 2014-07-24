% [U D V] = svd(L);
null_sol = V(:,end-1);
%     singular_values = diag(D);
    p=zeros(yn,xn,zn);
    vy=zeros(yn,xn,zn);
    vx=zeros(yn,xn,zn);
    vz=zeros(yn,xn,zn);
    
    % Process all Grid points
    for i=1:1:yn
        for j=1:1:xn
            for k = 1:1:zn
                % Global index for P, vx, vy in S()
                inp=((k-1)*xn*yn + (j-1)*yn+i)*4-3; % P
                invx=inp+1;
                invy=inp+2;
                invz=inp+3;
                % P
                p(i,j,k)=null_sol(inp)*kcont;
                % vx
                vx(i,j,k)=null_sol(invx);
                % vy
                vy(i,j,k)=null_sol(invy);
%                 vz
                vz(i,j,k)=null_sol(invz);
            end
        end
    end
    

[ax ay az] = meshgrid(1:xn,1:yn,1:zn);
vec3DToVtk(vx,vy,vz,ax,ay,az,'nullSol.vtk');



% To examine min and max values of v and p for specified range of modes.
% E.g for all the null space basis.
% start_indx = 5228;
% end_indx = 5324;
% vx_max = (1e-20)*ones(1,end_indx-start_indx);
% vy_max = vx_max;
% vz_max = vx_max;
% p_max = vx_max;
% vx_min = ones(1,end_indx-start_indx);
% vy_min = vx_min;
% vz_min = vx_min;
% p_min = vx_min;
% for i_indx = 1:end_indx-start_indx
%     null_sol = V(:,i_indx-1+start_indx);
% %     singular_values = diag(D);
%     p=zeros(yn,xn,zn);
%     vy=zeros(yn,xn,zn);
%     vx=zeros(yn,xn,zn);
%     vz=zeros(yn,xn,zn);
%     
%     % Process all Grid points
%     for i=1:1:yn
%         for j=1:1:xn
%             for k = 1:1:zn
%                 % Global index for P, vx, vy in S()
%                 inp=((k-1)*xn*yn + (j-1)*yn+i)*4-3; % P
%                 invx=inp+1;
%                 invy=inp+2;
%                 invz=inp+3;
%                 % P
%                 p(i,j,k)=null_sol(inp)*kcont;
%                 % vx
%                 vx(i,j,k)=null_sol(invx);
%                 % vy
%                 vy(i,j,k)=null_sol(invy);
% %                 vz
%                 vz(i,j,k)=null_sol(invz);
%             end
%         end
%     end
%     
%     vx_max(i_indx) = max(abs([vx(:); vx_max(i_indx)]));
%     vy_max(i_indx) = max(abs([vy(:); vy_max(i_indx)]));
%     vz_max(i_indx) = max(abs([vz(:); vz_max(i_indx)]));
%     p_max(i_indx) = max(abs([p(:); p_max(i_indx)]));
%     vx_min(i_indx) = min(abs([vx(:); vx_min(i_indx)]));
%     vy_min(i_indx) = min(abs([vy(:); vy_min(i_indx)]));
%     vz_min(i_indx) = min(abs([vz(:); vz_min(i_indx)]));
%     p_min(i_indx) = min(abs([p(:); p_min(i_indx)]));
% end
