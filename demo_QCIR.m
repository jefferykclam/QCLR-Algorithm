% This script demonstrates how the QCIR algorithm works.
% The function QCIR registers between the moving and static images with
% landmark constraints.
% The outputs of the QCIR algorithm are the deformation homeomorphism, the
% corresponding Beltrami coefficient, and the registered image.
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
% Example 1
% load('QCIR_Circleexample.mat');
% Example 2
% load('QCIR_ICexample2.mat');
% Example 3
load('QCIR_Bigcircle_example.mat');

%% QCIR algorithm
tic
[map,map_mu,register] = QCIR(moving,static,landmark,target,'plot',1,'level',5);
toc