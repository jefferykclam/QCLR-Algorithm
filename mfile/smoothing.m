function smooth_mu = smoothing(mu,Smooth_Operator,Operator)

% This function smooths the magnitude of Beltrami coefficients
%
% Inputs
% mu : m x 1 Beltrami coefficients
% Smooth_Operator : the smoothing operator
% Operator : sets of operators
%
% Outputs
% smooth_mu : m x 1 smoothed Beltrami coefficients
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

vmu = Operator.f2v*mu;
nvmu = Smooth_Operator\abs(vmu);
vmu = nvmu.*(cos(angle(vmu))+sqrt(-1)*sin(angle(vmu)));
smooth_mu = Operator.v2f*vmu;

end