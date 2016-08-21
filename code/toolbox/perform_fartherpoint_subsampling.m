function I = perform_fartherpoint_subsampling(D,N,k0, options)

% perform_fartherpoint_subsampling - perform subsampling by FP sampling
%
%   I = perform_fartherpoint_subsampling(D,N,k0);
%
%   D is a pairwise distance matrix.
%   I is a list of N indexes.
%   k0 is the initial point (optional).
%
%   Copyright (c) 2015 Gabriel Peyre

options.null = 0; 
verb = getoptions(options, 'verb', 0);

if nargin<3
    k0 = 1;
end

I = [k0];
d = min(D(k0,:),[],1);
for i=1:N-1
    if verb
        progressbar(i,N-1);
    end
    [~,i] = max( d ); i = i(1); 
    d = min(d, min(D(i,:),[],1) );
    I(end+1) = i;
end
I = I(:);

end