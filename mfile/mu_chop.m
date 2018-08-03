function mu_new = mu_chop(mu,bound,constant)

% This function hanges the norm of the input mu to be less than some constant
%
% Inputs
% mu : m x 1 Beltrami coefficients
% bound : upper limit of the norm of the Beltrami coefficient.
% constant : constant replace the norm of the Beltrami coefficients which exceed the upper limit
%
% Outputs
% mu_new : m x 1 updated Beltrami coefficient
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

mu_new = mu;
ind = abs(mu)>=bound;
mu_new(ind) = constant*(cos(angle(mu(ind)))+1i*sin(angle(mu(ind))));

end