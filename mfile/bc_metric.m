function mu = bc_metric(face, vertex, map, dimension)

% This function bc_metric computes the Beltrami coefficients of the
% mapping from "vertex" to "map", i.e. f(vertex) = map.
% Two choices of dimension can be chosen, which are 2 and 3 respectively.
% This dimension indicates whether the target map is planar or a 3D
% surface.
%
% Inputs
% face : m x 3 triangulation connectivity
% vertex : n x 3 vertices coordinates
% mu : m x 1 Beltrami coefficients
% dimension : 2 (planar mesh) or 3 (3D surface).
%
% Outputs
% A : The corresponding Beltrami coefficients.
%
% Function is written by Jeffery Ka Chun Lam (2013)
% www.jefferykclam.com
% Reference : 
% L. M. Lui, K. C. Lam, T. W. Wong and X. Gu, 
% Texture map and video compression using Beltrami representation.
% SIAM Journal on Imaging Sciences, 6(4):1880-1902, 2013.

face = face';
vertex = vertex';
map = map';


MyDifferential = diff_operator(vertex, face);

switch dimension
    case 2
        f = complex(map(1,:), map(2,:));
        mu = (MyDifferential.Dc*(f.')) ./ (MyDifferential.Dz*(f.'));
    case 3
        dXdu = transpose(MyDifferential.Dx* map(1,:)');
        dXdv = transpose(MyDifferential.Dy* map(1,:)');
        
        dYdu = transpose(MyDifferential.Dx* map(2,:)');
        dYdv = transpose(MyDifferential.Dy* map(2,:)');
        
        dZdu = transpose(MyDifferential.Dx* map(3,:)');
        dZdv = transpose(MyDifferential.Dy* map(3,:)');
        
        E = dXdu.^2 + dYdu.^2 + dZdu.^2;
        G = dXdv.^2 + dYdv.^2 + dZdv.^2;
        F = dXdu.*dXdv + dYdu.*dYdv + dZdu.*dZdv;

        mu = (E - G + 2 * 1i * F) ./ (E + G + 2*sqrt(E.*G - F.^2));
        mu = mu.';
    otherwise
        error('Dimension must be 3 or 2.')
end


end

function MyDifferential = diff_operator(vertex, face)
    n = size(face, 2);
    Mi = reshape([1:n;1:n;1:n], [1,3*n]);
    Mj = reshape(face, [1,3*n]);

    [e1,e2,e3]= GetEdge(vertex,face);
    
    area = GetSignedArea_edge(e1, e2);
    area = [area;area;area];

    Mx = reshape([e1(2,:);e2(2,:);e3(2,:)]./area /2 , [1, 3*n]);
    My = -reshape([e1(1,:);e2(1,:);e3(1,:)]./area /2 , [1, 3*n]);

    Dx = sparse(Mi,Mj,Mx);
    Dy = sparse(Mi,Mj,My);
    Dz = (Dx - 1i*Dy) / 2; Dc = (Dx + 1i*Dy) / 2;

    MyDifferential = struct('Dx', Dx, 'Dy', Dy, 'Dz', Dz, 'Dc', Dc);

end

function [e1, e2, e3] = GetEdge(vertex,face)
    e1 = vertex(1:2,face(3,:)) - vertex(1:2,face(2,:));
    e2 = vertex(1:2,face(1,:)) - vertex(1:2,face(3,:));
    e3 = vertex(1:2,face(2,:)) - vertex(1:2,face(1,:));
end

function area = GetSignedArea_edge(e1, e2)
    xb = -e1(1,:);
    yb = -e1(2,:);
    xa = e2(1,:);
    ya = e2(2,:);
    area = ( xa.*yb - xb.*ya )/2;
end
