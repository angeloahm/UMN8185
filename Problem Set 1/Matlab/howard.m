V0 = zeros(n_k,1);
V = zeros(n_k,1);
V_temp = zeros(n_k,1);
policy_k = zeros(n_k,1);
policy_c = zeros(n_k,1);
policy_index = zeros(n_k,1);

for iter=1:max_iter
    
    for i_k = 1:n_k
        
        k = k_grid(i_k);
        for i_k_prime = 1:n_k
            k_prime = k_grid(i_k_prime);
            c = A*k^alpha + (1-delta)*k - k_prime;
            if c<0
                c=1e-100;
            end
            V_temp(i_k_prime,1) = log(c) + beta*V0(i_k_prime,1);
        end
        [value, i_policy] = max(V_temp);
        V(i_k) = value;
        policy_index(i_k) = i_policy;
        policy_k(i_k) = k_grid(i_policy);
        policy_c(i_k) = A*k^alpha + (1-delta)*k - k_grid(i_policy);
    end
    
    for fake_iter = 1:m
        V_howard = log(policy_c) + beta*V(policy_index);
        V = V_howard;
    end
    
    
    norm = max(abs(V(:)-V0(:)));
    if norm<eps
        fprintf('VFI converged at iteration: %.5f \n', iter)
        break
    else
        fprintf('Error = %.8f; Iter = %.8f \n', norm, iter)
        V0 = V;
    end
    
end
