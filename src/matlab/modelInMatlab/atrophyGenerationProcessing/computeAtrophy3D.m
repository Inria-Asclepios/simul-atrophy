%  computeAtrophy
function a = computeAtrophy3D(m,n,r)

a = zeros(m,n,r);
% a = 2*ones(m,n,r);
% a(round(m/2),round(n/2),round(r/2)) = 0.2;


% a(round(m/2+m/4),round(n/2+n/4),round(r/2+r/4)) = 2;  %delta function
% a(round(m/2+m/4),round(n/2+n/4), round(n/2-n/4)) = 2;  %delta function
% a(round(m/2+m/4),round(n/2-n/4), round(n/2+n/4)) = 2;  %delta function
% a(round(m/2+m/4),round(n/2-n/4), round(n/2-n/4)) = 2;  %delta function
% a(round(m/2-m/4),round(n/2+n/4), round(n/2+n/4)) = 2;  %delta function
% a(round(m/2-m/4),round(n/2+n/4), round(n/2-n/4)) = 2;  %delta function
% a(round(m/2-m/4),round(n/2-n/4), round(n/2+n/4)) = 2;  %delta function
% a(round(m/2-m/4),round(n/2-n/4), round(n/2-n/4)) = 2;  %delta function


sw = 1;
a(round(m/2)-sw:round(m/2)+sw,round(n/2)-sw:round(n/2)+sw, round(r/2)-sw:round(r/2)+sw) = 0.2;  %small step function

% range = @(m ,sw) round(m/2+m/4)-sw:round(m/2+m/4)+sw; 
% a(range())

% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw,round(r/2+r/4)-sw:round(r/2+r/4)+sw) = 2; 


% gw = 40;
% a = 2*gbell([m n],[gw m-gw gw n-gw],[gw+20 gw+20], 2); %gaussian
% b = 2*gbell([m n],[gw m-gw gw n-gw],[gw-20 gw-20], 2); %gaussian
% a = a+b;
%
% a = (2*(exp([0:1/(100-1):1])))';
% a = -ceil((n-1)/2):floor((n-1)/2);
% sigma = 5;
% a = (exp(-(a.*a)/sigma))';

