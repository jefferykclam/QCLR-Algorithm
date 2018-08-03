function [map,map_mu] = QCLR(face,vertex,landmark,target,dimension,varargin)

% 	QCLR Algorithm
%
% Input : 
% face : m x 3 triangulation connectivity
% vertex : n x 3 vertices coordinates
% landmark : vertices ID as landmarks in the form of k x 1 vector
% target : corresponding target position in 2D in the form of k x 2 matrix
% dimension : dimension of the target rectangular domain in the form of [Height, Width]
% varargin
% Corner : corner vertexi indexes in anti-clockwise order
% mu : User defined initial Beltrami coefficient in the form of mx1 vector
% Parameter : Parameters input for the algorithm
% Plot : enabling to plot the result
%
% Output :
% map : registration result
% map_mu : corresponding Beltrami coefficients
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

cornerflag = 0;
parameterflag = 0;
plotflag = 0;
mu = 0*face(:,1);
IterationNumber = 0;

for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if ~isempty(regexp(varargin{i},'^[Cc]orner','match'))
        corner = varargin{i+1};
        cornerflag = 1;
    end
    if ~isempty(regexp(varargin{i},'^[Pp]arameter','match'))
        P = varargin{i+1};
        parameterflag = 1;
    end
    if ~isempty(regexp(varargin{i},'^[Mm]u','match'))
        mu = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Pp]lot','match'))
        plotflag = 1;
    end
end

if size(landmark,2) == 2
    landmark = vertex_search(landmark,vertex);
end

if size(landmark,1) ~= size(target,1)
	error('Please check landmark & target input!');
end

if ~cornerflag
    corner = vertex_search([...
        min(vertex(:,1)),min(vertex(:,2));
        max(vertex(:,1)),min(vertex(:,2));
        max(vertex(:,1)),max(vertex(:,2));
        min(vertex(:,1)),max(vertex(:,2))],vertex);
end

if ~parameterflag
    load('parameter_QCLR.mat');
end

Operator = createOperator(face,vertex);

[map,update_mu] = lbs_rect(face,vertex,mu,'Landmark',landmark,...
    'Target',target,'Corner',corner,'Height',dimension(1),'Width',dimension(2));
count = overlap(update_mu,P.upperBound);
tolerance = max(abs(update_mu - mu));
landmark_error = max(sum(abs(vertex(landmark,1:2)-target),2));

while (tolerance > .1) || ...
        (landmark_error > P.landmarkTolerance) || count ~= 0
    
    mu = update_mu;
    penalty = P.sigma;
    penalty = penalty + P.sigmaIncrease;
    
    Smooth_Operator = (1/(penalty))*(P.alpha*speye(length(vertex)) + ...
        penalty*speye(length(vertex)) - Operator.laplacian/2);
    smooth_mu = smoothing(update_mu,Smooth_Operator,Operator);
    smooth_mu = mu_chop(smooth_mu,P.upperBound,P.chopVal);
    
    for p = 1:P.dt_Number
        [~,update_mu] = lbs_rect(face,vertex,smooth_mu,'Landmark',landmark,...
            'Target',target,'Corner',corner,'Height',dimension(1),'Width',dimension(2));
        smooth_mu = update_mu + P.dt*(smooth_mu - update_mu);
        smooth_mu = mu_chop(smooth_mu,P.upperBound,P.chopVal);
    end
    
    [map,update_mu] = lbs_rect(face,vertex,smooth_mu,'Corner',corner,...
        'Height',dimension(1),'Width',dimension(2));
    
    landmark_error = max(sum(abs(map(landmark,1:2)-target),2));
    tolerance = max(abs(update_mu - mu));
    count = overlap(update_mu,P.upperBound);
    
    IterationNumber = IterationNumber + 1;
      
    if IterationNumber > P.iterationNumber
        warning('Iteration exceeds 100. Automatically terminated');
        break;
    end

end

map_mu = update_mu;

if plotflag
    % moving image
    show_mesh(face,vertex); 
    hold on; title('Moving');
    scatter3(vertex(landmark,1),vertex(landmark,2),vertex(landmark,3), 'ro', 'filled');
    scatter3(target(:,1),target(:,2),1+0*target(:,1), 'bo', 'filled');

    % registration result
    show_mesh(face,map); 
    title('Target');
end

end
    
    

