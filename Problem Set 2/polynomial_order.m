function p = polynomial_order(max_p, toler, alpha, c_l, c_u, utility, c_points)
    
for p = 1:max_p
    
   c_fine = linspace(c_l, c_u, c_points);

   if utility==1
       u_real = log(c_fine);
   else
       if utility==2
           u_real = c_fine.^(1/2);
       else
           u_real = c_fine.^(1-alpha) / (1-alpha);
       end
   end
   
   
   [coef, S] = polyfit(c_fine, u_real, p); 
   cheb_u = polyval(coef, c_fine);
    
   
   
   norm = max(abs((u_real(:) - cheb_u(:)) ./ u_real(:)));
   
   if norm<toler
       fprintf('\n')
       fprintf('Polynomial of order:   %.1f \n',p)
       break
   else
       fprintf('Norm   %.4f \n',norm)
       
   end
    
end




