% This script demonstrates how the QCLR algorithm works.
% The function QCLR registers two shapes by using landmark constraints.
% The outputs of the QCLR algorithm are the deformation homeomorphism and the
% corresponding Beltrami coefficien.
%
% PLEASE NOTE THAT THIS SCRIPT IS A SIMPLIFIED VERSION
% A FULL VERSION MAY BE RELEASED IN THE FUTURE
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.


addpath(genpath('example'));
addpath(genpath('mfile'));

%% loading example data
load('QCLR_example.mat');

%% QCLR algorithm
dimension = [50,50];
[map,map_mu] = QCLR(face,vertex,landmark,target,dimension,'plot',1);
