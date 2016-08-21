function [e,deriv,deriv_fd] = check_gradient(func, dims, nrepl, epsilon, tol)

% check_gradient - check gradient by comparison with finite differences
%
%   [e,deriv,deriv_fd] = check_gradient(func, dims, nrepl);
%
%   Copyright (c) 2015 Gabriel Peyre


if nargin<3
    nrepl = 10; 
end
if nargin<4
    epsilon = 1e-9;
end
if nargin<5
    tol = Inf;
end

if length(dims)==1
    dims = [dims 1];
end

dotp = @(x,y)sum(x(:).*y(:));

deriv = []; deriv_fd = [];
for i=1:nrepl
    x = randn(dims);
    d = randn(dims);
    [f,g] = func(x+epsilon*d);
    % true derivative
    deriv(i) = dotp( g,d );
    % finite differences
    [f1,~] = func(x);
    deriv_fd(i) = ( f-f1 ) / epsilon;
end
e = abs(deriv-deriv_fd)./max(abs(deriv));

if max(e)>tol
    warning('Problem with gradient');
end    
end