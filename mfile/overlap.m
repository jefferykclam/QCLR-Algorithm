function count = overlap(mu,bound)

% This function counts the number of |mu| greater than 1 or the input bound
%
% Inputs
% mu : m x 1 Beltrami coefficients
% bound : upper limit of the norm of the Beltrami coefficient.
%
% Outputs
% count : number of |mu| greater than 1 or the input bound
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

if nargin < 2
    count = sum(abs(mu) >= 1);
else
    count = sum(abs(mu) >= bound);
end

end