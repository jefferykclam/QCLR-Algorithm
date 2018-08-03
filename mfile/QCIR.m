function [map,map_mu,register,ms_map,ms_mu,ms_register] = QCIR(moving,static,landmark,target,varargin)

% QCIR Algorithm
%
% Inputs 
% moving : the moving image
% statuc : the static image
% landmark : vertices ID as landmarks in the form of k x 1 vector
% target : corresponding target position in 2D in the form of k x 2 matrix
% varargin
% Level : multiscale level (Default : 3)
% Plot : enabling to plot the result
%
% Output
% map : registration result
% map_mu : corresponding Beltrami coefficients
% register : registered image
%
% Function is written by Jeffery Ka Chun Lam (2013)
% www.jefferykclam.com
% Reference : 
% L. M. Lui, K. C. Lam, T. W. Wong and X. Gu, 
% Texture map and video compression using Beltrami representation.
% SIAM Journal on Imaging Sciences, 6(4):1880-1902, 2013.

level = 3;
plotflag = 0;

for i = 1:length(varargin)
    if isa(varargin{i},'double') 
        continue;
    end
    if ~isempty(regexp(varargin{i},'^[Ll]evel','match'))
        level = varargin{i+1};
    end
    if ~isempty(regexp(varargin{i},'^[Pp]lot','match'))
        plotflag = 1;
    end
end

[face,vertex] = image_meshgen(size(static,1),size(static,2));

if size(landmark,2) == 2
    landmark = vertex_search(landmark,vertex);
end

if size(landmark,1) ~= size(target,1)
	error('Please check landmark & target input!');
end

if (sum(size(moving)-size(static))~=0)
    moving = imresize(moving,size(static),'bicubic');
end

ms_moving = cell(level,1);
ms_static = cell(level,1);
ms_landmark = cell(level,1);
ms_target = cell(level,1);
ms_face = cell(level,1);
ms_vertex = cell(level,1);
ms_mu = cell(level,1);
ms_map = cell(level,1);
ms_register = cell(level,1);

for k = level:-1:1
    
    display(['Registering in Level ', int2str(k),'...']);
    
    ms_moving{k} = imresize(moving,...
        [round(size(static,1)/(2^(k-1))),round(size(static,2)/(2^(k-1)))],...
        'bicubic');
    ms_static{k} = imresize(static,...
        [round(size(static,1)/(2^(k-1))),round(size(static,2)/(2^(k-1)))],...
        'bicubic');
    [tempface, tempvertex] = image_meshgen(size(ms_static{k},1),size(ms_static{k},2));
    ms_face{k} = tempface;
    ms_vertex{k} = tempvertex;
    [temp_landmark,temp_index] = unique(vertex_search(vertex(landmark,:)/(2^(k-1)),tempvertex));
    ms_landmark{k} = temp_landmark;
    temp_target = target/(2^(k-1));
    ms_target{k} = temp_target(temp_index,:);
end

[initial_map,initial_mu] = QCLR(ms_face{level},ms_vertex{level},ms_landmark{level},...
    ms_target{level},size(ms_static{level}));

targetsize = size(ms_static{level});
initial_register = warpimage(ms_moving{level},targetsize,ms_vertex{level},initial_map);

for k = level:-1:2
    
    gaussian_size = [max(round(size(ms_static{k},1)/3),60),max(round(size(ms_static{k},2)/3),60),max(round(size(ms_static{k},1)/12),10)];
    
    [temp_map,temp_mu,temp_register] = registration_lm(ms_face{k},ms_vertex{k},initial_map,...
        initial_register,ms_static{k},ms_landmark{k},ms_target{k},gaussian_size);

    [temp_map,temp_mu] = QCLR(ms_face{k},ms_vertex{k},ms_landmark{k},...
    ms_target{k},size(ms_static{k}),'mu',temp_mu);

    targetsize = size(ms_static{k});
    temp_register = warpimage(ms_moving{k},targetsize,ms_vertex{k},temp_map);
    
    ms_register{k} = temp_register;
    ms_map{k} = temp_map;
    ms_mu{k} = temp_mu;
    
    interpolant_X = scatteredInterpolant(ms_vertex{k}(:,1),ms_vertex{k}(:,2),temp_map(:,1));
    interpolant_Y = scatteredInterpolant(ms_vertex{k}(:,1),ms_vertex{k}(:,2),temp_map(:,2));
    
    temp_vertex = ms_vertex{k-1};
    temp_vertex(:,1) = temp_vertex(:,1) * (max(ms_vertex{k}(:,1)) - min(ms_vertex{k}(:,1))) ...
        / (max(temp_vertex(:,1)) - min(temp_vertex(:,1))) * ( 1 - 1e-10);
    temp_vertex(:,2) = temp_vertex(:,2) * (max(ms_vertex{k}(:,2)) - min(ms_vertex{k}(:,2))) ...
        / (max(temp_vertex(:,2)) - min(temp_vertex(:,2))) * ( 1 - 1e-10);
    
    clear initial_map;
    initial_map(:,1) = interpolant_X(temp_vertex(:,1),temp_vertex(:,2));
    initial_map(:,2) = interpolant_Y(temp_vertex(:,1),temp_vertex(:,2));
    
    initial_map(:,1) = initial_map(:,1) * (max(ms_vertex{k-1}(:,1)) - min(ms_vertex{k-1}(:,1))) ...
        / (max(initial_map(:,1)) - min(initial_map(:,1))) * ( 1 - 1e-10);
    initial_map(:,2) = initial_map(:,2) * (max(ms_vertex{k-1}(:,2)) - min(ms_vertex{k-1}(:,2))) ...
        / (max(initial_map(:,2)) - min(initial_map(:,2))) * ( 1 - 1e-10);
    
    targetsize = size(ms_moving{k-1});
    initial_register = warpimage(ms_moving{k-1},targetsize,ms_vertex{k-1},initial_map);
    drawnow;
    
end
    
[map,map_mu,register] = registration_lm(ms_face{1},ms_vertex{1},initial_map,...
        initial_register,ms_static{1},ms_landmark{1},ms_target{1},gaussian_size);
    
[map,map_mu] = QCLR(ms_face{1},ms_vertex{1},ms_landmark{1},...
    ms_target{1},size(ms_static{1}),'mu',map_mu);

targetsize = size(ms_static{1});
register = warpimage(ms_moving{1},targetsize,ms_vertex{1},map);

ms_map{1} = map;
ms_register{1} = register;
ms_mu{1} = map_mu;

if plotflag
    figure
    imshow(moving); title('Moving Image with landmarks');
    hold on;
    if size(landmark,2) == 1
        plotlandmark = vertex(landmark,1:2);
        scatter3(plotlandmark(:,1),plotlandmark(:,2), 1+0*plotlandmark(:,1), 'ro', 'filled');
    else
        scatter3(landmark(:,1),landmark(:,2), 1+0*landmark(:,1), 'ro', 'filled');
    end

    figure
    imshow(static); title('Static Image with landmarks;');
    hold on;
    scatter3(target(:,1),target(:,2),1+0*target(:,1), 'bo', 'filled');

    figure
    imshow(register); title('Registration result');
end

end
    
    

