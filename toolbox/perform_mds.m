function [X,s] = perform_mds(D,d, method, options)

% perform_classical_mds - perform Multidimensional Scaling.
%
%   [X,s] = perform_mds(D,d, method, options);
%
%   method is either 'classical' or 'smacof'
%   D is a pairwise (N,N) distance matrix.
%   X is (N,d) embedding matrix.
%   s is the stress value during iterations (only for SMACOF)
%
%   Copyright (c) 2015 Gabriel Peyre

options.null = 0;
if nargin<3
    method = 'classical';
end

N = size(D,1);

switch method
    case 'classical'
        N = size(D,1);
        % centering
        M = -.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2));
        % spectral embedding
        opt.disp = 0; opt.isreal = 1; opt.issym = 1;
        [X, val] = eigs(M, d, 'LR', opt);
        % scaling
        for i=1:d
            X(:,i) = X(:,i)*sqrt(val(i,i));
        end   
         s = [];
    case 'smacof'
        X0 = getoptions(options, 'X0', randn(N,d));
        niter = getoptions(options, 'niter', 20);   
        verb = getoptions(options, 'verb', 0);
        Stress = @(d)sqrt( sum( abs(D(:)-d(:)).^2 ) / N^2 ); 
        DMat = @(Y)sqrt( repmat(sum(Y.^2),N,1) + repmat(sum(Y.^2),N,1)' - 2*Y'*Y);
        % do a global rescale
        Y = X0';
        D0 = DMat(Y);
        Y = Y * mean(D(:))/mean(D0(:)); 
        % operator B
        remove_diag = @(b)b - diag(sum(b));
        B = @(D1)remove_diag( -D./max(D1,1e-10) );
        % SMACOF algorithm
        Y = Y-repmat(mean(Y,2), [1 N]);
        s = [];
        for i=1:niter
            if verb
                progressbar(i,niter);
            end
            Y = Y * B(DMat(Y))' / N;
            % update
            Y = Y-repmat(mean(Y,2), [1 N]);
            % record stress
            s(end+1) = Stress(DMat(Y));
        end
        X = Y';        
    otherwise
        error('Unknown method');
end

X = real(X);

end