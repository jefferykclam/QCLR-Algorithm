function [M,Tx,Ty] = simple_demon(I1,I2,times,gaussian_size)

% This function computes the vector field using additive demon's algorithm.
%
% Inputs
% I1 : moving image
% I2 : static image
% times : number of moves of additive demons
% gaussian_size : size of the guassian kernel
%
% Outputs
% M : updated moving image
% Tx : x-coordinate vector field
% Ty : y-coordinate vector field
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

S=I2; M=I1;
alpha=2.5;

Hsmooth=fspecial('gaussian',[gaussian_size(1) gaussian_size(2)],gaussian_size(3)); % [60,60],10
Tx=zeros(size(M)); Ty=zeros(size(M));
[Sy,Sx] = gradient(S);
for itt=1:times
        Idiff=M-S;

        [My,Mx] = gradient(M);
        Ux = -Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        Uy = -Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2))+(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
 
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;
        Tx=Tx+3*imfilter(Ux,Hsmooth);
        Ty=Ty+3*imfilter(Uy,Hsmooth);

        D(:,:,1) = Ty;
        D(:,:,2) = Tx;
        M = imwarp(I1,D);
end