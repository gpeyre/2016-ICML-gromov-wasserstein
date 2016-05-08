function [Dbary,R] = compute_gw_barycenters(D,D0,lambda,epsilon,options)

% compute_gw_barycenters - compute GW barycenters of distance matrices
%
%   [Dbary,R] = compute_gw_barycenters(D,lambda,epsilon,options);
%
%   R display the evolution of the minimized error.
%
%   Copyright (c) 2016 Gabriel Peyre

options.null = 0;
bary_method = getoptions(options, 'bary_method', 'alternating');
gw_loss = getoptions(options, 'gw_loss', 'l2');
N0 = size(D{1},1);
normalize = @(x)x/sum(x(:));

%%
% Use unform densities
% for barycenter
mu = normalize(ones(size(D0,1),1));
% for the other
Mu = {};
for k=1:length(D)
    Mu{k} = normalize(ones(size(D{k},1),1));
end

switch bary_method
    case 'bfgs'
        %%% L-BFGS -- depreciated %%
        opts.verb = getoptions(options, 'bfgs_verb', 1);
        opts.niter = getoptions(options, 'bfgs_iter', 20);
        opts.bfgs_memory = getoptions(options, 'bfgs_memory', 10);
        opts.report = @(f,v)v;
        %
        func = @(d)LossCallback_Barycenters(D,Mu,d,mu, lambda, epsilon, options);
        func1 = @(x)VectorizeCallback(func,x,size(D{1}));
        [d, R, info] = perform_bfgs(func1, D0(:), opts);
        d = reshape(d,size(D{1}));
        % set diagonal to 0
        Dbary = d-diag(diag(d));
    case {'alternating' 'alterating'}
        %%% Alternating minimization %%%
        Dbary =D0; 
        tol_alternating = getoptions(options, 'tol_alternating', 1e-4);
        niter_alternating = getoptions(options, 'niter_alternating', 15);
        verb_alternating = getoptions(options, 'verb_alternating', 1);
        R = [];
        gamma = {};
        for i=1:niter_alternating
            if verb_alternating
                progressbar(i,niter_alternating);
            end
            Dnew = zeros(size(Dbary));
            R(i) = 0;
            for s=1:length(D)
                if i==1
                    options.gamma_init = Mu{s}*mu';
                else
                    options.gamma_init = gamma{s};
                end
                options.gamma_init = Mu{s}*mu'; %% Using clever initialiation is not a good idea apparently
                [gamma{s},a,b,GW] = perform_gw_sinkhorn(D{s}, Dbary, Mu{s}, mu, epsilon, options);
                % refine
                options.gamma_init = gamma{s}; 
                [gamma{s},a,b,GW] = perform_gw_sinkhorn(D{s}, Dbary, Mu{s}, mu, 0, options);                
                switch gw_loss
                    case 'l2'
                        Dnew = Dnew + lambda(s) * gamma{s}' * D{s} * gamma{s};
                    case 'kl'
                        Dnew = Dnew + lambda(s) * gamma{s}' * log( 1e-10 + D{s} ) * gamma{s};
                    otherwise
                        error('Unknown loss.');
                end
                R(i) = R(i) + lambda(s) * GW(end);
            end
            Dbary = Dnew ./ (mu*mu');
            if strcmp(gw_loss, 'kl')
                Dbary = exp(Dbary);
            end
            % early stopping
            if (i>2) && (tol_alternating>0) && ( (min(R(1:end-1)) - R(end))/abs(R(1))<tol_alternating )
                % stop if not enough progress
                if verb_alternating==1
                    progressbar(niter_alternating,niter_alternating);
                end
                break;
            end
        end
    otherwise
        error('Unknown');
end

end