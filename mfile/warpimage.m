function warp = warpimage(image, targetsize, initial_mesh, deformed_mesh)

% This function warps the moving image by the deformation determined by the
% original mesh and the deformed mesh.
% Inputs
% image : moving image
% targetsize : the target warping image size
% initia_mesh : the initial mesh embedding to the moving image
% deformed_mesh : the deformed mesh which describes the deformation
% Outputs
% warp : the warped image
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.
    

inverse_vector = initial_mesh(:,1:2) - deformed_mesh(:,1:2);
Fx = scatteredInterpolant(deformed_mesh(:,1),deformed_mesh(:,2),inverse_vector(:,1));
Fy = scatteredInterpolant(deformed_mesh(:,1),deformed_mesh(:,2),inverse_vector(:,2));
ivector(:,1) = Fx(initial_mesh(:,1),initial_mesh(:,2));
ivector(:,2) = Fy(initial_mesh(:,1),initial_mesh(:,2));
D(:,:,1) = (reshape(ivector(:,1),targetsize));
D(:,:,2)  = (reshape(ivector(:,2),targetsize));
warp = imwarp(image,D,'linear');

end