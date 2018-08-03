function [map,map_mu,Edge] = lbs_rect(face,vertex,mu,varargin)

% This function lbs restores the quasi-conformal map in rectangular domain 
% corresponding to the given Beltrami coefficients mu by minimizing the least 
% square energy from the Beltrami equation.
% Corners corresponding to the target rectangular domain can be chosen
% manually.
% Hard landmark (both interior / Dirichlet boundary) constraints can be
% added, where landmarks are the vertices index, target position is represented 
% by the 2-dimensional coordinates. 
%
% Inputs
% face : m x 3 triangulation connectivity
% vertex : n x 3 vertices coordinates
% mu : m x 1 Beltrami coefficients
% varargin :
% Corner : corner vertexi indexes in anti-clockwise order
% Landmark : vertices ID as landmarks in the form of k x 1 vector
% Target : corresponding target position in 2D in the form of k x 2 matrix
%
% Outputs
% map : corresponding quasi-conformal map f
% map_mu : The true, admissible Beltrami coefficients
% Edge : Divided edges
%
% Remarks on Corner input
%
% Defualt :
% 	4 extreme points are chosen.
%
% User input :
% 	Anti-clockwise location for the corners.
%   	P4 ---(Edge 3)--- P3
%   	|	       		   |
%	(Edge 4)			(Edge 2)
%   	|	               |
%   	P1 ---(Edge 1)--- P2
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

heightflag = 0;
cornerflag = 0;
width = 1;
landmark = [];
target = [];

for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if ~isempty(regexp(varargin{i},'^[Cc]orner','match'))
        corner = varargin{i+1};
        cornerflag = 1;
    end
    if ~isempty(regexp(varargin{i},'^[Hh]eight','match'))
        height = varargin{i+1};
        heightflag = 1;
    end
    if ~isempty(regexp(varargin{i},'^[Ww]idth','match'))
        width = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Ll]andmark','match'))
        landmark = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Tt]arget','match'))
        target = varargin{i+1};
    end
    
end

if size(landmark,1) ~= size(target,1)
	error('Please check landmark & target input!');
end

if isempty(landmark)
    targetlmX = [];
    targetlmY = [];
else
    targetlmX = target(:,1);
    targetlmY = target(:,2);
end
    
[Ax,abc,area] = generalized_laplacian2D(face,vertex,mu); Ay = Ax;
bx = zeros(length(vertex),1); by = bx;
if ~cornerflag
    corner = vertex_search([...
        min(vertex(:,1)),min(vertex(:,2));
        max(vertex(:,1)),min(vertex(:,2));
        max(vertex(:,1)),max(vertex(:,2));
        min(vertex(:,1)),max(vertex(:,2))],vertex);
end
Edge = close_curve_division(freeBoundary(triangulation(face,vertex)),corner);

VBdyC = [Edge{4};Edge{2}]; VBdy = [Edge{4}*0; Edge{2}*0 + width];
landmarkx = [landmark;VBdyC]; targetx = [targetlmX;VBdy];
bx(landmarkx) = targetx;
Ax(landmarkx,:) = 0; 
Ax(landmarkx,landmarkx) = diag(ones(length(landmarkx),1));
mapx = Ax\bx;

if ~heightflag
    f = face';
    F = f(:);
    ux = mapx(F(1:3:end-2)).*(vertex(F(2:3:end-1),2) - vertex(F(3:3:end),2)) -...
        mapx(F(2:3:end-1)).*(vertex(F(1:3:end-2),2) - vertex(F(3:3:end),2)) +...
        mapx(F(3:3:end)).*(vertex(F(1:3:end-2),2) - vertex(F(2:3:end-1),2));
    uy = mapx(F(2:3:end-1)).*(vertex(F(1:3:end-2),1) - vertex(F(3:3:end),1)) -...
        mapx(F(1:3:end-2)).*(vertex(F(2:3:end-1),1) - vertex(F(3:3:end),1)) -...
        mapx(F(3:3:end)).*(vertex(F(1:3:end-2),1) - vertex(F(2:3:end-1),1));
    height = 0.25*sum((abc(:,1).*ux.^2 + 2*abc(:,2).*ux.*uy + abc(:,3).*uy.^2)./area);
end

HBdyC = [Edge{1};Edge{3}]; HBdy = [Edge{1}*0; Edge{3}*0 + height];
landmarky = [landmark;HBdyC]; targety = [targetlmY;HBdy];
by(landmarky) = targety;
Ay(landmarky,:) = 0; 
Ay(landmarky,landmarky) = diag(ones(length(landmarky),1));
mapy = Ay\by;

map = [mapx,mapy,0*vertex(:,1)+1];
map_mu = bc_metric(face,vertex,map,2);

end