function [gamma,a,b,GW] = perform_gw_sinkhorn(d1,d2,mu1,mu2,epsilon, options)

% perform_gw_sinkhorn - perform iterative sinkhorn iterations to solve for GW matching
%
%   [gamma,a,b,GW] = perform_gw_sinkhorn(d1,d2,mu1,mu2,epsilon, options);
%
%   Solves for (a stationary point of)
%       GW = min_{gamma} C - <gamma*d2,d1*gamma> + epsilon*Entropy(gamma)
%       s.t.   gamma*1=mu1 and gamma^T*1=mu2  
%   where C = 1/2*(<d1*mu1,mu1> + <d2*mu2,mu2>) is a constant offset
%
%   Set options.log_domain=1 to use log-domain Sinkhorn (stable even with small epsilon). 
%
%   Copyright (c) 2015 Gabriel Peyre

% GW = @(gamma,d1,d2)-dotp(d1*gamma,gamma*d2);

options.null = 0;
niter = getoptions(options, 'niter', 50);
niter_sinkhorn = getoptions(options, 'niter_sinkhorn', 400);
display_mode = getoptions(options, 'display_mode', 0); 
tol_gw = getoptions(options, 'tol_gw', 0); 
verb = getoptions(options, 'verb', 1);
log_domain = getoptions(options, 'log_domain', 1);
safe = getoptions(options, 'safe', 0);

gw_loss = getoptions(options, 'gw_loss', 'l2');

dotp = @(x,y)sum(x(:).*y(:));
Entropy = @(x)sum(x(:).*(log(x(:)+1e-20)-1) );

% constant offset of the cost
switch gw_loss
    case 'l2'
        f1 = @(a)a.^2/2;
        f2 = @(b)b.^2/2;
        h1 = @(a)a;
        h2 = @(b)b;
    case 'kl'            
        f1 = @(a)a.*log(a+1e-15) - a;
        f2 = @(b)b;
        h1 = @(a)a;
        h2 = @(b)log(b+1e-15);
    otherwise 
        error('Unknown loss.');
end

GWcst = dotp( f1(d1)*mu1,mu1) + dotp( f2(d2)*mu2,mu2);

GW = [];
% initial match
gamma = getoptions(options, 'gamma_init', []);
if isempty(gamma)
    gamma = mu1*mu2';
end

opts.niter = niter_sinkhorn;
opts.tol = getoptions(options, 'tol_sinkhorn', 0);
opts.niter_sink = getoptions(options, 'niter_log_sink', 1); % log-exp switch frequency
opts.verb = 0;
opts.safe = safe;


for i=1:niter
    if verb
        progressbar(i,niter);
    end
	C = - h1(d1)*gamma*h2(d2');
    if 0
        C0 = C;
    else
        C0 = C-min(C(:));
    end
    % run sinkhorn
    if log_domain==0
        [gamma,a,b] = perform_sinkhorn(C0,epsilon, mu1,mu2, opts);
        opts.a_init = a; % warm restart
    else
        % run sinkhorn, log domain
        [gamma,r,s] = perform_sinkhorn_log(C0,epsilon, mu1,mu2,opts);
        opts.r_init = r; opts.s_init = s;
    end
    %
    if 0
        % plot dual mutipliers
        clf; hold on;
        plot(log(a)-mean(log(a)), 'r');
        plot(log(b)-mean(log(b)), 'b');
    end
    GW(end+1) = GWcst + dotp( C,gamma ) + epsilon * Entropy(gamma);
    if (i>2) && (tol_gw>0) && ( (min(GW(1:end-1)) - GW(end))/abs(GW(1))<tol_gw )
        % stop of not enough progress
        if verb==1
            progressbar(niter,niter);
        end
        break;
    end
    if display_mode
        clf;
        imagesc(gamma);
        drawnow;
    end
    if isnan(gamma(1)),
        error('Computation blew up, epsilon too small');
    end
end

if not(exist('a'))
    a = exp(r/epsilon); b = exp(s/epsilon);
end

end