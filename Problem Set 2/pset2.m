%% Problem Set 2
% Angelo Mendes

clear 
%% Comparing interpolation methods
alpha = 2;  %CRRA Parameter
theta = 1;  %Equally spaced grid
n = 10;     %Gridpoints
c_u = 5;    %Upper limit for the grid of c
c_l = 0.1;  %Lower limit for the grid of c

%Create naive grid
c = linspace(c_l, c_u, n);

%Polinomially expanded 
c_expanding = c_l + (c_u - c_l) .* c .^ theta;
% Finer grid
c_fine = linspace(min(c_expanding), max(c_expanding), 10000);

%Store values for the utility functions
u_1 = log(c_expanding);
u_2 = c_expanding.^(1/2);
u_3 = c_expanding.^(1-alpha) / (1-alpha);

%Store values for the utility functions
u1_real = log(c_fine);
u2_real = c_fine.^(1/2);
u3_real = c_fine.^(1-alpha) / (1-alpha);


u1_linear = interp1(c_expanding, u_1, c_fine);
u1_spline = interp1(c_expanding, u_1, c_fine, 'spline');

u2_linear = interp1(c_expanding, u_2, c_fine);
u2_spline = interp1(c_expanding, u_2, c_fine, 'spline');

u3_linear = interp1(c_expanding, u_3, c_fine);
u3_spline = interp1(c_expanding, u_3, c_fine, 'spline');

h = figure;
plot(c_fine, u1_real, c_fine, u1_linear, '--', c_fine, u1_spline, '--')
xlabel('c')
ylabel('u(c)')
title('Interpolation methods - log(c)')
legend('Function', 'Linear', 'spline')
saveas(h,sprintf('u1_interp_theta%d.png',theta))

h = figure;
plot(c_fine, u2_real, c_fine, u2_linear, '--', c_fine, u2_spline, '--')
xlabel('c')
ylabel('u(c)')
title('Interpolation methods - c^{1/2}')
legend('Function', 'Linear', 'spline')
saveas(h,sprintf('u2_interp_theta%d.png',theta))

h = figure;
plot(c_fine, u3_real, c_fine, u3_linear, '--', c_fine, u3_spline, '--')
xlabel('c')
ylabel('u(c)')
title(sprintf('Interpolation methods - CRRA (alpha =%d)',alpha))
legend('Function', 'Linear', 'spline')
saveas(h,sprintf('u3_interp_theta%d.png',theta))

%% Using function to find optimal grindpoints
n_guess = 2;
toler = 0.01;
[n, norm] = optimal_gridpoints(n_guess, toler, 10, 4, c_l, c_u, 'linear', 3);

%% Extrapolation 
% Take the coefficients for for the interpolations and evaluate off gridpoints
alpha=5;
extrapolation_grid = [0.05 5.1 5.5 6];

extrapolation_u1_linear = interp1(c, log(c), extrapolation_grid, 'linear', 'extrap');
extrapolation_u2_linear = interp1(c, c.^(1/2), extrapolation_grid, 'linear', 'extrap');
extrapolation_u3_linear = interp1(c, c.^(1-alpha) ./ (1-alpha), extrapolation_grid,...
    'linear', 'extrap');

extrapolation_u1_spline = interp1(c, log(c), extrapolation_grid, 'spline', 'extrap');
extrapolation_u2_spline = interp1(c, c.^(1/2), extrapolation_grid, 'spline', 'extrap');
extrapolation_u3_spline = interp1(c, c.^(1-alpha) ./ (1-alpha), extrapolation_grid,...
    'spline', 'extrap');


% An alternative way of doing that would be the following:
%
% extrapolation_u1_spline = zeros(4,1);
% extrapolation_u2_spline = zeros(4,1);
% extrapolation_u3_spline = zeros(4,1);
%
% Save spline coefs:
% spline_coef_u1 = spline(c, log(c));
% spline_coef_u2 = spline(c, c.^(1/2));
% spline_coef_u3 = spline(c, c.^(1-alpha) ./ (1-alpha));
% 
% Evaluate off-gridpoints using the spline coefs: 
% for i=1:4
%    c_extra = extrapolation_grid(i); 
%    
%    extrapolation_u1_spline(i) = ppval(spline_coef_u1, c_extra);
%    extrapolation_u2_spline(i) = ppval(spline_coef_u2, c_extra);
%    extrapolation_u3_spline(i) = ppval(spline_coef_u3, c_extra);
%      
% end

%% Chebyshev Polynomials

% p1 = poly(c);
% y1 = polyval(p1, c_fine);
p = 5; %Degree of the Chebyshev polynomial
%Coefs of the polynomial that better fits (in a LS sense) data in u1_real
[coef1, S1] = polyfit(c_fine, u1_real, p); 
cheb_u1 = polyval(coef1, c_fine);

[coef2, S2] = polyfit(c_fine, u2_real, p); 
cheb_u2 = polyval(coef2, c_fine);



h = figure;
plot(c_fine, u1_real, c_fine, cheb_u1)
xlabel('c')
ylabel('u(c)')
title('Chebyshev Polynomial Approximation')
legend('Function', 'Polynomial')
saveas(h, 'u1_cheb', 'png')

h = figure;
plot(c_fine, u2_real, c_fine, cheb_u2)
xlabel('c')
ylabel('u(c)')
title('Chebyshev Polynomial Approximation')
legend('Function', 'Polynomial')
saveas(h, 'u2_cheb', 'png')

alpha = 10;  %CRRA Parameter
u3_real = c_fine.^(1-alpha) / (1-alpha);
[coef3, S3] = polyfit(c_fine, u3_real, p); 
cheb_u3 = polyval(coef3, c_fine);

h = figure;
plot(c_fine, u3_real, c_fine, cheb_u3)
xlabel('c')
ylabel('u(c)')
title('Chebyshev Polynomial Approximation')
legend('Function', 'Polynomial')
saveas(h,sprintf('u3_cheb_alpha%d.png',alpha))
%% Optimal polynomial order
p = polynomial_order(500, 0.01, 2, c_l, c_u, 2, 1000);











