function [E,GradE] = VectorizeCallback(func,X,sz)

% VectorizeCallback - helper to use an arbitrary shaped callback into e.g. L-BFGS
%
%   [E,GradE] = VectorizeCallback(func,X,sz);
%
%   func takes input of size sz.
%   GradE is of size (prod(sz),1), i.e. it is a column vector.
%
%   Copyright (c) 2015 Gabriel Peyre

X = reshape(X,sz);
[E,GradE] = func(X);
GradE = GradE(:);

end