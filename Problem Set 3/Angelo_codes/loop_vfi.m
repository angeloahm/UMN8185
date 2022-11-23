V0 = zeros(n_A, n_k);
V = zeros(n_A, n_k);
V_temp = zeros(n_k,1);
c_policy = zeros(n_A, n_k);
k_policy = zeros(n_A, n_k);
policy_index = zeros(n_A, n_k);

utility = @(c) (c .^ (1-gamma)) / (1-gamma);

for i_A=1:n_A
    for i_k=1:n_k
        A = A_grid(i_A);
        k = k_grid(i_k);
        
        c = A * (k^alpha) + (1-delta)*k;
                
        if c>0
            u = (c^(1-gamma)) / (1-gamma);
        else
            u = -1e+7;
        end
        V0(i_A, i_k) = u;
    end
end

for iter=1:max_iter
    
    for i_A=1:n_A
        A = A_grid(i_A);
        for i_k=1:n_k
            k = k_grid(i_k);
            for i_k_prime=1:n_k
                k_prime = k_grid(i_k_prime);
                
                % Compute consumption
                c = A * (k^alpha) + (1-delta)*k - k_prime;
                
                if c>0
                    u = (c^(1-gamma)) / (1-gamma);
                else
                    u = -1e+7;
                end
                
                V_temp(i_k_prime) = u + beta*dot(P(i_A,:), V0(:,i_k_prime));
                
            end
            [value, ind] = max(V_temp);
            V(i_A, i_k) = value;
            policy_index(i_A, i_k) = ind;
            k_policy(i_A, i_k) = k_grid(ind);
            c_policy(i_A, i_k) = A*k^alpha + (1-delta)*k - k_grid(ind);
            
        end
    end
    
%     for fake_iter = 1:m
%        V_howard = utility(c_policy) + beta*(P*V(policy_index));
%        V = V_howard;
%     end
        
    norm = max(abs(V(:)-V0(:)));
    
    if norm<eps
        fprintf('VFI converged at iteration: %.5f \n', iter)
        break
    else
        fprintf('Error = %.8f; Iter = %.8f \n', norm, iter)
        V0 = V;
    end
    
end
    
    