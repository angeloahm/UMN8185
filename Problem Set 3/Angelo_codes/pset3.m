%% Problem Set 3
% Angelo Mendes

% Preliminaries

clear

% Define function
g = @(p) 0.5*p^(-0.5) + 0.5*p^(-0.2);
% Real derivative
derivative = @(p) -0.5*0.5*p^(-1.5)-0.2*0.5*p^(-1.2);
eps_grid = 10 .^ (-linspace(1,10,10)); %Tolerance grid

%% Exercise 1 - One Sided Derivative

% Define the numerical derivative
g_prime = @(p,eps) (g(p+eps)-g(p)) / eps;


g_prime_grid = zeros(length(eps_grid), 1);

for i=1:length(g_prime_grid)
    
    eps = eps_grid(i);
    g_prime_grid(i) = g_prime(1.5, eps);
    
end

error_1s = g_prime_grid - derivative(1.5);
[error, i_eps_optimal] = min(abs(error_1s(:)));
fprintf('Minimum difference is: %d\n', error)
fprintf('Best epsilon is: %d\n', eps_grid(i_eps_optimal))


%% Exercise 2 - Two Sided Derivative

g_prime_2s = @(p, eps) (g(p+eps) - g(p-eps))/(2*eps);
g_prime_2s_grid = zeros(length(eps_grid), 1);

for i=1:length(g_prime_grid)
    
    eps = eps_grid(i);
    g_prime_2s_grid(i) = g_prime_2s(1.5, eps);
    
end

error_2s = g_prime_2s_grid - derivative(1.5);
[error, i_eps_optimal] = min(abs(error_2s(:)));
fprintf('Minimum difference is: %d\n', error)
fprintf('Best epsilon is: %d\n', eps_grid(i_eps_optimal))

%% Exercise 3 - Root finding
f = @(p) g(p)-0.75;
p0_bisection = bisection(f, 1, 1355, 100, 1e-6);
p0_secant = secant(f, 1, 1355, 100, 1e-6);
p0_newton = newton(f, derivative, 1, 100, 1e-6);
p0_fzero = fzero(f, [1,1355]);

%% Exercise 4 - Brent's method

%% Exercise 5

clear; clc;
% AR(1) parameters
sigma = 0.007;                                              % SD of the A shock
rho = 0.98;                                             % Persistence of the A shock
n_A = 15;                                               % A gridpoints
[log_A_grid, P] = rouwenhorst(n_A, 0, rho, sigma);      % Discretize A grid
A_grid = exp(log_A_grid);

% Numerical parameters
theta = 1.5;                                            % Spacing the grid
n_k = 71;                                               % Capital gridpoints 
eps_derivative = 1e-5;                                  % Numerical derivative tolerance
max_iter = 5000;                                        % Max iter (VFI)
m = 30;                                                  % Howard iterations
eps = 1e-5;                                             % Tolerance VFI

                  
% Model parameters
alpha=0.36;                                              % Production parameter
gamma=2;                                                % CRRA parameter
delta=0.7;                                              % Depreciation
beta=0.95;                                              % Discount factor


k_upper = (alpha*max(A_grid(:)) ./ (1/beta + delta -1)) .^ (1/(1-alpha));
k_lower = (alpha*min(A_grid(:)) ./ (1/beta + delta -1)) .^ (1/(1-alpha))*0.1;
kss = (alpha / (1/beta + delta -1)) .^ (1/(1-alpha));

% Capital grid
%k_grid = k_lower + (k_upper - k_lower) .* linspace(k_lower, k_upper, n_k) .^ theta;   
k_grid = linspace(0.75*kss,1.2*kss,n_k);

tic
vfi
toc

k_policy_int = k_policy;

tic
loop_vfi
toc

k_policy_loop = k_policy;

h = figure;
plot(k_grid, k_policy_loop(1,:), k_grid, k_policy_loop(15,:))
xlabel('k')
ylabel('k^{\prime}(A,k)')
title('Policy Functions')
legend()

hold on 
plot(k_grid, k_policy_int(1,:), '--', k_grid, k_policy_int(15,:), '--')
legend('A^{VFI}_{min}', 'A^{VFI}_{max}', 'A^{Int}_{min}', 'A^{Int}_{max}')





