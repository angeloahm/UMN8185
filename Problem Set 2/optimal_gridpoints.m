function [n, norm] = optimal_gridpoints(n_guess, toler, alpha, theta, c_l, c_u, method, utility)
    
for iter = 1:10000
    
   n=n_guess;
   c = linspace(c_l, c_u, n);
   c_expanding = c_l + (c_u - c_l) .* c .^ theta;
   c_fine = linspace(min(c_expanding), max(c_expanding), 10000);

   if utility==1
       u = log(c_expanding);
       u_real = log(c_fine);
   else
       if utility==2
           u = c_expanding.^(1/2);
           u_real = c_fine.^(1/2);
       else
           u = (c_expanding.^(1-alpha)) / (1-alpha);
           u_real = c_fine.^(1-alpha) / (1-alpha);
       end
   end
   
   
   u_linear = interp1(c_expanding, u, c_fine);
   u_spline = interp1(c_expanding, u, c_fine, 'spline');
    
   
   if method=='linear'
      u_inter = u_linear;
   else
      u_inter = u_spline;
   end
   
   norm = max(abs((u_real(:) - u_inter(:)) ./ u_real(:)));
   
   if norm<toler
       fprintf('\n')
       fprintf('Optimal n is:   %.1f \n',n)
       break
   else
       fprintf('Norm   %.4f \n',norm)
       
   end
   n_guess = n+1;
    
end




