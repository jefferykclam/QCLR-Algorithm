function index = vertex_search(XYZ,vertex)

% This function searches vertex indexes for which the coordinates of the
% corresponding vertices are nearest to the input XY / XYZ.
%
% Inputs
% XYZ : (X,Y)/(X,Y,Z) coordinates of ponits in the form of k x 2 / k x 3
% matrices
% vertex : n x 3 vertices coordinates
%
% Output :
%	index : k x 1 vertex indexes of the points (X,Y)/(X,Y,Z) 
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

if size(XYZ,2) ~= 2 && size(XYZ,2) ~= 3
    error('Input feature points must be a kx2 or kx3 matrix');
end


k = size(XYZ,1);
n = size(vertex,1);
index = zeros(k,1);
v = vertex';

switch size(XYZ,2)
    case 2
		[~,index] = min(sqrt((repmat(v(1,:),k,1)-repmat(XYZ(:,1),1,n)).^2 +...
            (repmat(v(2,:),k,1)-repmat(XYZ(:,2),1,n)).^2),[],2);
    case 3
		[~,index] = min(sqrt((repmat(v(1,:),k,1)-repmat(XYZ(:,1),1,n)).^2 +...
            (repmat(v(2,:),k,1)-repmat(XYZ(:,2),1,n)).^2 +...
            (repmat(v(3,:),k,1)-repmat(XYZ(:,3),1,n)).^2),[],2);
end

if length(unique(index)) ~= length(index)
    warning('Some points are sharing the same vertex index found')
end   

end