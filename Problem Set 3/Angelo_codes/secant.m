function root = secant(f, x0, x1, max_iter, toler)

for iter = 1:max_iter
    x_new = x1 - ((x1-x0)/(f(x1)-f(x0)))*f(x1);
    
    if abs(f(x_new))<toler
        root=x_new;
        fprintf('Secant concluded in: %d iterations \n', iter)
        fprintf('Root is:                           %d \n', root)
        fprintf('Error:                             %d \n', abs(f(root)))
        break
    else
        fprintf('Error:                             %d \n', abs(f(x_new)))
        x1=x_new;
    end
end
end