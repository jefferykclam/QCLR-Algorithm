function [map,map_mu] = lbs(face,vertex,mu,landmark,target)

% This function lbs restores the quasi-conformal map corresponding to
% the given Beltrami coefficients mu by minimizing the least square energy
% from the Beltrami equation
% Hard landmark (both interior / Dirichlet boundary) constraints can be
% added, where landmarks are the vertices index, target position is represented 
% by the 2-dimensional coordinates. 
%
% Inputs
% face : m x 3 triangulation connectivity
% vertex : n x 3 vertices coordinates
% mu : m x 1 Beltrami coefficients
% landmark : vertices ID as landmarks in the form of k x 1 vector
% target : corresponding target position in 2D in the form of k x 2 matrix
%
% Outputs
% map : corresponding quasi-conformal map f
% map_mu : The true, admissible Beltrami coefficients
%
% Function is written by Jeffery Ka Chun Lam (2013)
% www.jefferykclam.com
% Reference : 
% L. M. Lui, K. C. Lam, T. W. Wong and X. Gu, 
% Texture map and video compression using Beltrami representation.
% SIAM Journal on Imaging Sciences, 6(4):1880-1902, 2013.

targetc = target(:,1) + 1i*target(:,2);
A = generalized_laplacian2D(face,vertex,mu);
b = -A(:,landmark)*targetc;
moving = setdiff(linspace(1,size(vertex,1),size(vertex,1)),landmark);
map(landmark,1) = targetc;
map(moving,1) = A(moving,moving)\b(moving);

map = [real(map),imag(map),1+0*map];
map_mu = bc_metric(face,vertex,map,2);

end