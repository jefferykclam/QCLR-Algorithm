function [face,vertex] = image_meshgen(height,width)

% This function image_meshgen creates mesh with given image size
%
% Input
% height : number of pixels in Y-axis (height)
% width : number of pixels in X-axis (width)
%
% Outputs
% face : m x 3 trangulation connectivity
% vertex : n x 3 vertices coordinates
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.


[x,y] = meshgrid(0:height-1, 0:width-1);
vertex = [x(:),y(:),1+0*x(:)];
face = delaunay(x(:),y(:));

end