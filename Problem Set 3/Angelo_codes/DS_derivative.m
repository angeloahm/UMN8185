function derivative = DS_derivative(f, x, eps)

derivative = (f(x+eps) - f(x-eps)) / (2*eps);


end