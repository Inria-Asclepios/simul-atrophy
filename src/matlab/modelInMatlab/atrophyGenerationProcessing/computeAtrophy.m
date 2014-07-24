% computeAtrophy
function a  = computeAtrophy(m,n)

a = zeros(m,n);
% a = 0.2*ones(m,n);
% a(ceil(m/2),ceil(n/2)) = 0.4;
a(ceil(m/2)+1,ceil(n/2)+1) = 0.2;
a(ceil(m/2)+1,ceil(n/2)-1) = 0.2;
a(ceil(m/2)-1,ceil(n/2)+1) = 0.2;
a(ceil(m/2)-1,ceil(n/2)-1) = 0.2;

% a(ceil(m/2),ceil(n/2)-2) = 0.2;
% a(ceil(m/2),ceil(n/2)+2) = 0.2;
% a(ceil(m/2)-2,ceil(n/2)) = 0.2;
% a(ceil(m/2)+2,ceil(n/2)) = 0.2;

% a(round(m/2+m/4),round(n/2+n/4)) = -0.2;  %delta function
% a(round(m/2+m/4),round(n/2-n/4)) = -0.2;  %delta function
% a(round(m/2-m/4),round(n/2-n/4)) = -0.2;  %delta function
% a(round(m/2-m/4),round(n/2+n/4)) = -0.2;  %delta function

% sw = 1;
% a(round(m/2)-sw:round(m/2)+sw,round(n/2)-sw:round(n/2)+sw) = 0.2;  %small step function

% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw) = 2; 
% a(round(m/2+m/4)-sw:round(m/2+m/4)+sw,round(n/2-n/4)-sw:round(n/2-n/4)+sw) = 2; 
% a(round(m/2-m/4)-sw:round(m/2-m/4)+sw,round(n/2-n/4)-sw:round(n/2-n/4)+sw) = 2; 
% a(round(m/2-m/4)-sw:round(m/2-m/4)+sw,round(n/2+n/4)-sw:round(n/2+n/4)+sw) = 2; 



% Centered Gaussian, 
% sigma = 1;
% a_scale = 0.9;
% [y x] = meshgrid(-floor(m/2):ceil(m/2)-1, -floor(n/2):ceil(n/2)-1);
% a = a_scale * exp(-1*(x.*x + y.*y)/sigma);
% a = ones(m,n);
% a(round(m/2+m/4),round(n/2-n/4)) = 2;  %delta function
% a(round(m/2-m/3),round(n/2-n/3)) = 2;  %delta function
% a(round(m/2-m/3):round(m/2+m/3),round(n/2-n/3):round(n/2+n/3)) = 2;   %step function

% gw = 10;
% a = 2*gbell([m n],[gw m-gw gw n-gw],[round(m/2) round(n/2)], 1); %gaussian


% gw = 20;
% a = 20*gbell([m n],[gw m-gw gw n-gw],[round(m/2) round(n/2)-5], 2); %gaussian
% b = 20*gbell([m n],[gw m-gw gw n-gw],[round(m/2) round(n/2)+5], 2); %gaussian
% a = a+b;

%
% a = (2*(exp([0:1/(100-1):1])))';
% a = -ceil((n-1)/2):floor((n-1)/2);
% sigma = 5;
% a = (exp(-(a.*a)/sigma))';


end

