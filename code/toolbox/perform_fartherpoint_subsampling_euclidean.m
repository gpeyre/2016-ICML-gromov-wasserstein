function I = perform_fartherpoint_subsampling_euclidean(X,N,k0)

% perform_fartherpoint_subsampling_euclidean - compute a equi-spaced subsampling
%
%   I = perform_fartherpoint_subsampling_euclidean(X,N,1) 
%
%   I is an index set of size N
%   X is of size (P,d)
%
%   Copyright (c) 2015 Gabriel Peyre

I = [k0];
d = distmat( X(I(end),:)', X' );
for i=1:N-1
    [~,I(end+1)] = max(d);
    d = min( distmat( X(I(end),:)', X' ), d );    
end
I = I(:);

end