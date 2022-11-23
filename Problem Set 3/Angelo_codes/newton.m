function root = newton(f, f_prime, x0, max_iter, toler)

for iter = 1:max_iter
    x_new = x0 - f(x0)/f_prime(x0);
    
    if abs(f(x_new))<toler
        root=x_new;
        fprintf('Newton-Raphson concluded in: %d iterations \n', iter)
        fprintf('Root is:                           %d \n', root)
        fprintf('Error:                             %d \n', abs(f(root)))
        break
    else
        fprintf('Error:                             %d \n', abs(f(x_new)))
        x0=x_new;
    end
end
end