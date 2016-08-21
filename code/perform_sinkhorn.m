function [gamma,a,b,Err] = perform_sinkhorn(C,epsilon,mu,nu, options)

% perform_sinkhorn - perform sinkhorn algorithm
%
%   [gamma,a,b,Err] = perform_sinkhorn(C,epsilon,mu,nu, options)
%
%   solve for (a,b) such that gamma = diag(a)*K*diag(b)
%   satisfies
%       gamma*ones(N,1) = mu
%       gamma'*ones(N,1) = nu
%
%   Err(i,:) is relative error on marginals at iterate i
%   options.niter controls #iterations
%   options.tol is for early stopping under prescribed accuracy on marginals
%
%   Copyright (c) 2015 Gabriel Peyre

if epsilon==0
    [~,gamma] = mexEMD(mu,nu,C);
    a = [];
    b = [];
    Err = [];
    return;
end

options.null = 0;
niter = getoptions(options, 'niter', 100);
tol = getoptions(options, 'tol', 1e-10);
a = getoptions(options, 'a_init', ones(size(C,1),1) );
safe = getoptions(options, 'safe', 0 );
% gibbs kernel
K = exp(-C/epsilon);

if safe
    K(K<1e-200)=1e-200; % Safe
end

Err = [];
for i=1:niter
    b = nu ./ (K'*a);
    if nargout>3 || tol>0
        Err(i,1) = norm(a.*(K*b)-mu, 1);
    end
    a = mu ./ (K*b);
    if nargout>3 || tol>0
        Err(i,2) = norm(b.*(K'*a)-nu, 1);
    end
    if max(Err(i,:))<tol
        break;
    end
end
gamma = diag(a)*K*diag(b);

end
