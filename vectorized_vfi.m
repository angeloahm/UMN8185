V0 = zeros(n_k,1);
V = zeros(n_k,1);
V_temp = zeros(n_k,1);
U = zeros(n_k, n_k);
policy_k = zeros(n_k,1);
policy_c = zeros(n_k,1);
policy_index = zeros(n_k,1);

for i_k=1:n_k
    for i_k_prime=1:n_k
        
        k = k_grid(i_k);
        k_prime = k_grid(i_k_prime);
        c = A*k^alpha + (1-delta)*k - k_prime;
        if c<0
            c=1e-10;
        end
        U(i_k, i_k_prime) = log(c);
        
    end
end


for iter=1:max_iter %start VFI

    V = U + beta*V0;






end
for iter=1:max_iter %start VFI
    
    for i_k = 1:n_k %for each k in the grid
        k = k_grid(i_k);
        j = 1*(i_k==1) + i_capital*(i_k>1);
        for i_k_prime = j:n_k %for each k_prime
            k_prime = k_grid(i_k_prime);
            %Define consumption
            c = A*k^alpha + (1-delta)*k - k_prime;
            if c<0
                c=1e-100;
            end
            
            % RHS of the Bellman equation
            V_temp(i_k_prime) = log(c) + beta*V0(i_k_prime);
            
            % We now use concavity and avoid typing "max" (check Cezar lecture notes)
            if i_k_prime<n_k
                if V_temp(i_k_prime)>V_temp(i_k_prime+1)
                    i_capital = i_k_prime; 
                    value = V_temp(i_k_prime);
                end
            else
                [value, i_capital] = max(V_temp);
            end
            
 
        end
        V(i_k) = value;
        policy_index(i_k) = i_capital;
        policy_k(i_k) = k_grid(i_capital);
        policy_c(i_k) = A*k^alpha + (1-delta)*k - k_grid(i_capital);
    end
    
    % Start policy iteration 
    for fake_iter = 1:m
        V_howard = log(policy_c) + beta*V(policy_index);
        V = V_howard;
    end
    
    %Define the Mac-Queen Porteus bounds (check Fatih's book for a good explanation)
    c_upper = (beta/(1-beta)) * max(V(:)-V0(:));
    c_lower = (beta/(1-beta)) * min(V(:)-V0(:));
    
    norm = c_upper - c_lower; %define the new norm
    if norm<eps
        V_median = V + ((c_upper+c_lower)/2);
        fprintf('VFI converged at iteration: %.5f \n', iter)
        break
    else
        fprintf('Error = %.8f; Iter = %.8f \n', norm, iter)
        V0 = V;
    end
    
end
