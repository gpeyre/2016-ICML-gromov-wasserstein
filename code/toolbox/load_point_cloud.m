function X = load_point_cloud(name, N, d)

% load_point_cloud - load a set of points 
%
%   X = load_point_cloud(name, N, d);
%
%   d is the dimension
%   N number of samples
%   X is of size (N,d)
%
%   Copyright (c) 2015 Gabriel Peyre

circle = @(N,c)[cos(2*pi*(0:N-1)'/N)+c(1),c(2)+ sin(2*pi*(0:N-1)'/N)];

switch name
    case 'gaussian'
        X = randn(N,d);
    case 'circle'
        X = circle(N, [0 0]);
    case 'grid'
        n = ceil(sqrt(N));
        t = linspace(-1,1,n)*.8;
        [b,a] = meshgrid(t,t);
        X = [a(:) b(:)];
    otherwise 
        error('Unknown cloud name');        
end

end