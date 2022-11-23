V0 = zeros(n_A, n_k);
V = zeros(n_A, n_k);
c_policy = zeros(n_A, n_k);
k_policy = zeros(n_A, n_k);
k_guess = zeros(n_A, n_k);

utility = @(c) (c .^ (1-gamma))/(1-gamma);

for i_A=1:n_A
    A = A_grid(i_A);
    for i_k=1:n_k
        k = k_grid(i_k);
        V0(i_A, i_k) = utility(A*k^alpha + (1-delta)*k);
        
    end
end

for iter=1:max_iter
    

    
    for i_A=1:n_A
        A = A_grid(i_A);
        
        inter = spline(k_grid, P(i_A,:)*V0);
        EV = @(kp) ppval(inter, kp);
        
        for i_k=1:n_k
            k = k_grid(i_k);
            %FOC = @(k_prime) (-1)*(A*k^alpha + (1-delta)*k - k_prime)^(-gamma)+ ...
            %    beta*DS_derivative(EV, k_prime, eps_derivative);
            VF = @(k_prime) -utility(A*k^alpha + (1-delta)*k - k_prime) - beta*EV(k_prime);
            ub = A*k^alpha + (1-delta)*k;
            lb = k_lower * 0.9;
            %k_prime = fzero(FOC, [lb ub]);
            options.Display = 'off';
            [k_prime, value] = fmincon(VF, k, [], [], [], [], lb, ub, [], options);
                
            k_policy(i_A, i_k) = k_prime;
            c_policy(i_A, i_k) = A*k^alpha + (1-delta)*k - k_prime;
            V(i_A, i_k) = value;
        end
    end 
    
    for fake_iter = 1:m
        V_howard = utility(c_policy) + beta*EV(k_prime);
        V = V_howard;
    end
    
    norm = max( (k_policy(:)-k_guess(:)) ./ (1+abs(k_policy(:))) );
    
    if norm<eps
        fprintf('VFI converged at iteration: %.5f \n', iter)
        break
    else
        fprintf('Error = %.8f; Iter = %.8f \n', norm, iter)
        V0 = V;
        k_guess = k_policy;
    end
    
    
    
end


        