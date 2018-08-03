function [map,map_mu,register] = registration_lm(face,vertex,previous_map,moving,static,landmark,target,gaussian_size)

% This is an internal function which updates the registration mapping by
% using both intensity and landmark constraints.
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

map = previous_map;
Operator = createOperator(face,vertex);
load('parameter_QCLR.mat');
temp_moving = moving;

if size(landmark,2) > 1
    landmark = vertex_search(landmark,vertex);
end

for k = 1:5

    [~,Bheight,Bwidth] = simple_demon(temp_moving,static,100,gaussian_size); 
    Iy = scatteredInterpolant(vertex(:,1),vertex(:,2),Bheight(:));
    Ix = scatteredInterpolant(vertex(:,1),vertex(:,2),Bwidth(:));
    
    Bheight = Iy(map(:,1),map(:,2));
    Bwidth = Ix(map(:,1),map(:,2));
    
    updated_map(:,1) = map(:,1) - Bwidth;
    updated_map(:,2) = map(:,2) - Bheight;
    
    mu = bc_metric(face,vertex,updated_map,2);
    
    P.alpha = 0.005;
    penalty = 0.5; % P.sigma
    penalty = penalty + P.sigmaIncrease;
    
    Smooth_Operator = (1/(penalty))*(P.alpha*speye(length(vertex)) + ...
        penalty*speye(length(vertex)) - Operator.laplacian/2);
    
    smooth_mu = smoothing(mu,Smooth_Operator,Operator);
    smooth_mu = mu_chop(smooth_mu,P.upperBound,P.chopVal);
    
    updated_map = lbs_rect(face,vertex,smooth_mu,'height',size(static,1)-1,'width',size(static,2)-1,'landmark',landmark,'target',target);
    mu = bc_metric(face,vertex,updated_map,2);
    
    inverse_vector = map(:,1:2) - updated_map(:,1:2);
    Fx = scatteredInterpolant(updated_map(:,1),updated_map(:,2),inverse_vector(:,1));
    Fy = scatteredInterpolant(updated_map(:,1),updated_map(:,2),inverse_vector(:,2));
    ivector(:,1) = Fx(vertex(:,1),vertex(:,2));
    ivector(:,2) = Fy(vertex(:,1),vertex(:,2));
    D(:,:,1) = (reshape(ivector(:,1),size(static)));
    D(:,:,2) = (reshape(ivector(:,2),size(static)));
    temp_moving = imwarp(temp_moving,D);

    map = updated_map;
    
end

map_mu = mu;
register = warpimage(moving, size(static), vertex, updated_map);

end
