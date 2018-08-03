function Edge = close_curve_division(B,pt)

% This function close_curve_division divides the close boundary into
% segments.
%
% Inputs
% B : n x 2 close Boundary edges 
%       For example B = [1 2; 2 3; 3 4; 4 1];
% pt : k x 1 points index for dividing the boundary
%
% Outputs
% Edge - MATLAB cell structure, segments of the close boundary
%
% Function is written by Jeffery Ka Chun Lam (2013)
% www.jefferykclam.com
% Reference : 
% L. M. Lui, K. C. Lam, T. W. Wong and X. Gu, 
% Texture map and video compression using Beltrami representation.
% SIAM Journal on Imaging Sciences, 6(4):1880-1902, 2013.

for i = 2:length(B)
    index1 = find(B(i:end,1)==B(i-1,2));
    if ~isempty(index1)
        tempa = B(i-1+index1,:);
        tempb = B(i,:);
        B(i,:) = tempa;
        B(i-1+index1,:) = tempb;
    else
        tempa = B(i-1+index1,:);
        tempb = B(i,:);
        B(i,:) = fliplr(tempa);
        B(i-1+index1,:) = tempb;
    end
end

n = length(pt);
Edge = cell(1,n);
[~,location] = ismember(pt,B(:,1));
[sort_location,index] = sort(location);

for i = 1:n-1
    Edge{i} = B(sort_location(i,1):sort_location(i+1,1),1);
end
Edge{n} = [B(sort_location(end,1):end,1);B(1:sort_location(1,1),1)];

Edge = Edge(index);

if min(cellfun(@length,Edge)) == 0
	error('Something wrong! Please check the input!')
end

end
