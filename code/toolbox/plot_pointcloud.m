function plot_pointcloud(X,col)

% plot_pointcloud - plot point cloud with colors
%
%   plot_pointcloud(X,col);
%
%   Copyright (c) 2015 Gabriel Peyre

if nargin<2
    col = X;
end

N = size(X,1);
ms = 10;
ms = 30;

for i=1:size(col,2)
    col(:,i) = rescale(col(:,i));
end
if size(col,2)==2
    col(:,3) = col(:,1)*0;
end

hold on;
for i=1:N
    if size(X,2)==3
        plot3(X(i,1), X(i,2), X(i,3), '.', 'MarkerSize', ms, 'color', col(i,:));
    else
        plot(X(i,1), X(i,2), '.', 'MarkerSize', ms, 'color', col(i,:));
    end
end
if size(X,2)==2
    axis ij;
end
axis tight; axis equal;  axis off; axis ij;

end