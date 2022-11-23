function [kgrid, T, VALUEFUNCTION, K_TOMORR, Kstar] = ngm(N)

%% Input

N;

%% Code

%%% Value Function Iteration (Stochastic and discrete grid)

%% Parameters
Par.alpha = 0.36;   % Capital participation
Par.gamma = 2;      % Coefficient of risk aversion
Par.beta  = 0.99;   % Discount factor
Par.delta = 0.03;   % Depreciation
Par.rho = 0.98;     % Persistence of the log of the productivity level
Par.sigma = 0.007;  % Standard deviation of shocks to the log of the productivity level
crit  =  1e-6;  %  Toleranz value

%% Solve for the steady state
Zstar = 0;
Kstar = ((1/Par.beta - 1 + Par.delta)./(Par.alpha * exp(Zstar))).^(1/(Par.alpha-1));

%% Create a grid for Z
meanZ = 0;
Grid.nZ = 15;  % number of points in our grid for Z
numStdZ = 2;  % number of standard deviations to cover with the grid
[Grid.Z, Grid.PZ]  = tauchen(Grid.nZ, meanZ, Par.rho, Par.sigma, numStdZ);

Grid.PZ = Grid.PZ'; % this is a 15 x 15 transition matrix for which the columns sum to 1
% the (i,j) element is the probability of moving from j to i.

%% Create a grid for K
%N = 31; %COMMENT FOR THE LOOP
Grid.nK = N; % Number of points in the K grid
Grid.K = linspace(0.75*Kstar, 1.25*Kstar,Grid.nK)';  % this is a 20 x 1 array of evenly spaced points

% %% Create a product of the two grids
% [ZZ,KK] =meshgrid(Grid.Z,Grid.K); %Create the Cartesian product of two vectors
% Grid.KK = KK(:);
% Grid.ZZ = ZZ(:);


%% Bellman iteration

kmin=min(Grid.K);
kmax=max(Grid.K);

% Out structures
kgrid=Grid.K;
zgrid=Grid.Z;
prob=Grid.PZ;
alpha=Par.alpha;
beta=Par.beta;
gamma=Par.gamma;
delta=Par.delta;

% Howard Improvement
H = 10;
%H = [1 5 10 25 100 500]; % Number of Howard steps

sizeH = size(H,2);

for howard=1:sizeH

tic

    h=H(howard);

% Loop Value Function Iteration

% Initial values
  val   = zeros(Grid.nK,Grid.nZ);
  tval  = zeros(Grid.nK,Grid.nZ);
  kdeci = zeros(Grid.nK,Grid.nZ);
  util  = zeros(Grid.nK,Grid.nK,Grid.nZ);
  gap =1; 
  n=1;
  
  while gap>=crit
      kdeci_old = kdeci;
    for i=1:Grid.nK
      for m=1:Grid.nZ
          k1 = fminbnd(@(k) valfun(k, kgrid, zgrid, alpha, gamma, beta, delta, i, m, val, prob),kmin,kmax); 
          tval(i,m) = -valfun(k1, kgrid, zgrid, alpha, gamma, beta, delta, i, m, val, prob);
          kdeci(i,m) = k1;
      end
    end

    if h>1
    val1 = tval;
    tval2=zeros(Grid.nK, Grid.nZ);
    for u=1:h  %howard steps
        for j=1:Grid.nK
            for z=1:Grid.nZ
            tval2(j,z) = -valfun(kdeci(j,z), kgrid, zgrid, alpha, gamma, beta, delta, j, z, val1, prob);
            end
        end
        val1=tval2;
    end
    %gap   = max(max(abs(tval-val)))           % 1 x 1  
    gap = max(max(abs((kdeci-kdeci_old)./(1+kdeci_old))));  % 1 x 1 
    val   = tval2;  % nK x nZ
    n=n+1;  
    
    else
    %gap   = max(max(abs(tval-val)))           % 1 x 1  
    gap = max(max(abs((kdeci-kdeci_old)./(1+kdeci_old))));  % 1 x 1 
    val   = tval;  % nK x nZ
    n=n+1; 
    end
  end

T(howard)=toc
      
% %% Results
% 
% DK = Grid.K/Kstar-1; % Capital grid as percent deviation from steady state
% 
% DKp = kdeci./Grid.K - 1;
% 
% %DKp = reshape(Kdeci,Grid.nK,Grid.nZ)./reshape(Grid.KK,Grid.nK,Grid.nZ) - 1;
%   % savings policy rule as a 20 x Z array expressed as a percent change from current K
% 
% v1=figure;
%     plot(DK, DKp);  % plot the policy rule
%     hold on;        % next plots go on the same figure
%     plot(DK, zeros(Grid.nK,1), 'k--'); % add a zero line, k-- means black and dashsed
%     xlabel('K in % deviation from steady state')  % label the axes
%     ylabel('(K'' - K)/K')
%     saveas(v1,sprintf('Howard_value_%d.png',h));
% 
% v2=figure; % Value Function
%     plot(Grid.K,val);
%     legend('Computational solution');
%     title('Value Function');
%     xlabel('kt'); 
%     ylabel('v(kt)');
% 
%     saveas(v2,sprintf('Howard_value2_%d.png',h));
% 
%  f2=figure; % Policy Function
%     plot(Grid.K,kdeci); hold;
%     plot(Grid.K, Grid.K, 'k');
%     legend('Computational solution','45 degree line');
%     title('Policy Function');
%     xlabel('kt'); 
%     ylabel('kt+1');
%     
%     saveas(f2,sprintf('Howard_policy%d.png',h));

end

%% Output
T;
VALUEFUNCTION = val(:,8); %% Keep the median shock
K_TOMORR = kdeci(:,8);
kgrid;
Kstar;

end