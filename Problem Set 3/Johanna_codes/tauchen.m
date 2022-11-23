function [Z,Zprob] = tauchen(nZ,a,rho,sigma,m)
%Function TAUCHEN
%
%Purpose:    Finds a Markov chain whose sample paths
%            approximate those of the AR(1) process
%                z(t+1) = a + rho * z(t) + eps(t+1)
%            where eps are normal with stddev sigma
%
%Format:     {Z, Zprob} = Tauchen(nZ,a,rho,sigma,m)
%
%Input:      nZ      scalar, number of points in z-grid
%            a       scalar
%            rho     scalar
%            sigma   scalar, std. dev. of epsilons
%            m       max +- std. devs.
%
%Output:     Z       1 * nZ vector, grid for Z
%            Zprob   nZ * nZ matrix, transition probabilities
%
%    Martin Floden
%    Fall 1996
%
%    This procedure is an implementation of George Tauchen's algorithm
%    described in Ec. Letters 20 (1986) 177-181.
%


Z     = zeros(nZ,1);
Zprob = zeros(nZ,nZ);

Z(nZ) = m * sqrt(sigma^2 / (1 - rho^2));
Z(1)  = -Z(nZ);
zstep = (Z(nZ) - Z(1)) / (nZ - 1);

for i=2:(nZ-1)
    Z(i) = Z(1) + zstep * (i - 1);
end 

Z = Z + a / (1-rho);

for j = 1:nZ
    for k = 1:nZ
        if k == 1
            Zprob(j,k) = cdf_normal((Z(1) - a - rho * Z(j) + zstep / 2) / sigma);
        elseif k == nZ
            Zprob(j,k) = 1 - cdf_normal((Z(nZ) - a - rho * Z(j) - zstep / 2) / sigma);
        else
            Zprob(j,k) = cdf_normal((Z(k) - a - rho * Z(j) + zstep / 2) / sigma) - ...
                         cdf_normal((Z(k) - a - rho * Z(j) - zstep / 2) / sigma);
        end
    end
end

Z = Z';


function c = cdf_normal(x)
    c = 0.5 * erfc(-x/sqrt(2));

