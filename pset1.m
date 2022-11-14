%% Problem Set 1 
% Angelo Mendes


%% Exercise 1 - Rounding Error
clear 

% Parameters
a=1;
b=100000;
n_points = 8;
gridn = -linspace(1,8,n_points);

roots = zeros(n_points,2);

for i=1:n_points
    
   n = gridn(i);
   c = 10^n;
   Delta = b^2 - 4*a*c;
   
   % First method
   x1 = (-b-(Delta)^(1/2)) / (2*a);
   x2 = (-b+(Delta)^(1/2)) / (2*a);
   roots(i, 1) = max(x1, x2);
   
   % Second  method
   q = (1/2)*(-b+sign(b)*(Delta)^(1/2));
   x3 = q/a;
   x4 = q/c;
   roots(i,2) = max(x3, x4);
   
    
end


plot(gridn, roots(:,1), gridn, roots(:,2));
xlabel('n')
ylabel('Larger root')
title('Larger root under different methods')
legend('Method 1', 'Method 2')

%% Exercise 2 - Truncation Error/Unstable Algorithms
clear 

phi = ((5)^(1/2)-1) / 2;
grid_phi = ones(20, 2);

grid_phi(2,2) = 0.61803398;

for n=2:21
    grid_phi(n,1) = (0.61803398)^(n-1);
end

for n=1:18
   grid_phi(n+2,2) = grid_phi(n,2) - grid_phi(n+1,2); 
    
end

plot(linspace(1,21,21), grid_phi(:, 1)-grid_phi(:, 2))

title('Difference between \phi under two different methods')

%xline(1) in more recent versions of MATLAB
xlabel('n')
ylabel('\phi')

%% Exercise 3 - Neoclassical Growth Model (VFI)

clear

% Parameters
A = 1; 
beta = 0.96;
delta = 1; %low delta ==> ~linear policies 
n_k = 1001;
alpha = 0.4;
eps = 1e-8;
max_iter = 5000;
m = 500;                    %Howard's parameter
fundamental_iterations = 2; %Good practice! Start with just VFI in the first couple of iterations

k_ss = (alpha*A / (1/beta + delta -1) ) ^ (1/(1-alpha));
k_grid = linspace(0.5*k_ss, 1.5*k_ss, n_k);


a = (1/(1-beta)) * (1/(1-alpha*beta)) * (log(A) + (1-alpha*beta)*log(1-alpha*beta)) ...
    + alpha*beta*log(alpha*beta);
b = (alpha / (1-alpha*beta));

k_analytical = (alpha*beta*A) * k_grid .^ alpha;
c_analytical = A*k_grid.^alpha - k_analytical;

tic
mp_bounds
toc

% if delta==1
%     figure(4)
%     plot(k_grid, k_analytical, k_grid, policy_k)
%     title('Comparison with Analytical Solution')
%     legend('Analytical','Numerical')
%     xlabel('k')
%     ylabel('g(k)')
%     
%     figure(5)
%     plot(k_grid, c_analytical, k_grid, policy_c)
%     title('Comparison with Analytical Solution')
%     legend('Analytical','Numerical')
%     xlabel('k')
%     ylabel('c(k)')
%     
% end
% 
% figure(1)
% plot(k_grid, V)
% title('Value function')
% xlabel('k')
% ylabel('V(k)')
% 
% figure(2)
% plot(k_grid, policy_k, k_grid, k_grid, '--')
% title('Capital Policy')
% legend('g(k)','45º line')
% xlabel('k')
% ylabel('g(k)')
% 
% figure(3)
% plot(k_grid, policy_c)
% title('Consumption Policy')
% xlabel('k')
% ylabel('c(k)')
%         
        







