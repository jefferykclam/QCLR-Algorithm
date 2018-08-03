function [map,map_mu] = QCLR_bdy(face,vertex,landmark,target,boundary,varargin)

% 	QCLR Algorithm
%
% QCLR algorithm specifies for rectangular domains
%
% Function is written by Jeffery Ka Chun Lam (2014)
% www.jefferykclam.com
% Reference : 
% K. C. Lam and L. M. Lui, 
% Landmark and intensity based registration with large deformations via Quasi-conformal maps.
% SIAM Journal on Imaging Sciences, 7(4):2364--2392, 2014.

parameterflag = 0;
plotflag = 0;
mu = 0*face(:,1);
IterationNumber = 0;

for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
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

if ~parameterflag
    load('parameter_QCLR.mat');
end

landmark = [landmark;boundary(:,1)];
target = [target;boundary(:,2:3)];

Operator = createOperator(face,vertex);

[map,update_mu] = lbs(face,vertex,mu,landmark,target);
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
        [~,update_mu] = lbs(face,vertex,smooth_mu,landmark,target);
        smooth_mu = update_mu + P.dt*(smooth_mu - update_mu);
        smooth_mu = mu_chop(smooth_mu,P.upperBound,P.chopVal);
    end
    
    [map,update_mu] = lbs(face,vertex,smooth_mu,landmark,target);
    
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
    show_mesh(face,vertex);title('Original');
    hold on
    scatter3(vertex(landmark,1),vertex(landmark,2),vertex(landmark,3),'ro','filled');
    show_mesh(face,map);title('Registered');
    hold on
    scatter3(target(:,1),target(:,2),map(3,landmark)','bo','filled');
    scatter3(map(landmark,1),map(landmark,2),map(landmark,3),'ro','filled');
end

end
    
    

