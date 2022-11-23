function tval=valfun(k, kgrid, zgrid, alpha, gamma, beta, delta, i, m, val, prob)

% This program gets the value function for a neoclassical growth model with
% uncertainty and CRRA utility

% do the interpolation
gg = interp1(kgrid,val,k,'cubic');

exp_gg = gg*prob(m,:)';

c = (exp(zgrid(m))*kgrid(i)^alpha + (1-delta)*kgrid(i) - k);  %  consumption(kt,kt+1, zt)

if c<0
    tval = -8888888888888888-800*abs(c);
else
    tval = (1/(1-gamma))*(c^(1-gamma)) + beta*exp_gg;
end
tval = -tval; % make it negative since we're maximizing and code is to minimize.