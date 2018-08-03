function show_mesh(face,vertex,coloring)

% This function show_mesh plots the mesh in either 2D / 3D surfaces.
%
% Inputs;
% face : m x 3 triangulation connectivity
% vertex : n x 3 vertices coordinates
% coloring : Color of each face / vertex. It can either be m x 3 matrices or
% n x 1 vectors.
%
% Function is written by Jeffery Ka Chun Lam (2013)
% www.jefferykclam.com
% Reference : 
% L. M. Lui, K. C. Lam, T. W. Wong and X. Gu, 
% Texture map and video compression using Beltrami representation.
% SIAM Journal on Imaging Sciences, 6(4):1880-1902, 2013.

if nargin < 3
    graph = figure; patch('Faces',face,'Vertices',vertex,'FaceColor','none',...
        'EdgeColor',[64,150,190]/255,'LineWidth',0.5);
else
    if ~isreal(coloring)
        coloring = abs(coloring);
    end
    graph = figure; patch('Faces',face,'Vertices',vertex,'FaceColor','flat','FaceVertexCData',coloring,...
        'EdgeColor','k', 'LineWidth',0.1);
    colorbar
end

axis equal
axis tight
axis vis3d
set(gcf,'Color','white');


assetData = struct('Vertex',vertex);
setappdata(gca,'AssetData',assetData);

dcm_obj = datacursormode(graph);
set(dcm_obj,'UpdateFcn',@ID_display,'Enable','on')
end

function txt = ID_display(obj,event_obj)

hAxes = get(get(event_obj,'Target'),'Parent');
assetData = getappdata(hAxes,'AssetData');

pos = get(event_obj,'Position');
id = vertex_search([pos(1),pos(2),pos(3)],assetData.Vertex);
txt = {['(x,y,z) : (',num2str(pos(1)),',',num2str(pos(2)),',',num2str(pos(3)),')'],...
    ['vertex ID : ',int2str(id)]};

end
