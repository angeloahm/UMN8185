[KK, AA, KK_prime] = meshgrid(k_grid, A_grid , k_grid);

c_policy = zeros(n_A, n_k);
k_policy = zeros(n_A, n_k);

utility = @(c) (c .^ (1-gamma))/(1-gamma);

CC = AA .* KK .^ alpha + (1-delta) .* KK - KK_prime;
U = utility(CC);
U(CC<0) = -1e7;

V0 = zeros(n_A, n_k);
V = zeros(n_A, n_k);
policy_index = zeros(n_A, n_k);

for iter=1:max_iter
    EV = repmat(P*V0, 1, 1, n_k);
    V_temp = U + beta*EV;
    [V, policy_index] = max(V_temp, [], 3);
    k_policy = k_grid(policy_index);

    norm = max(abs(V(:)-V0(:)));
    
    if norm<eps
        fprintf('VFI converged at iteration: %.5f \n', iter)
        break
    else
        fprintf('Error = %.8f; Iter = %.8f \n', norm, iter)
        V0 = V;
    end
    
end
    
    