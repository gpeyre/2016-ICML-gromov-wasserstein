function [gamma,r,s,err] = perform_sinkhorn_log(c,epsilon, p,q, options)

% perform_sinkhorn_log - perform Sinkhorn over the log-domain.
%
%   [gamma,r,s,err] = perform_sinkhorn_log(p,q, c,epsilon, options);
%
%   Perform Sinkhorn algorithm partially over the log-domain (i.e. with
%   stabilization) to solve
%       min_{gamma} <c,gamma> + epsilon*sum( gamma(:).*log(gamma(:)) )
%   subject to
%       gamma*1=p
%       gamma'*1=q
%
%   (r,s) are dual variables such that 
%       gamma = exp( -(c-r*1'-1*s')/epsilon );
%
%   err(i) = norm( gamma*1-p,1 ) for the current iterate i.
%
%   option.tol is a tolerance, forcing exit if err(i)<tol.
%
%   options.niter is the total number of iterations.
%   options.niter_sink is the number of inner Sinkhorn stabilization loop.
%      - Setting niter_sink=niter corresponds to the usual Sinkhorn, which is
%       unstable for small epsilon.
%      - Setting niter_sink=1 corresponds to doing all the iterates over
%       the log domain, which is slow.
%
%   Set options.verb=0 if you do not want to display a progressbar. 
%
%   Copyright (c) 2015 Gabriel Peyre


options.null = 0;
niter = getoptions(options, 'niter', 1000);
niter_sink = getoptions(options, 'niter_sink', 50);
tol = getoptions(options, 'tol', 0);
verb = getoptions(options, 'verb', 1);

[Np,Nq]=size(c);

up = ones(Np,1);  uq = ones(Nq,1); 

r = getoptions(options, 'r_init', zeros(Np,1) );  
s = getoptions(options, 's_init', zeros(Nq,1) );
err = [];
for i=1:ceil(niter/niter_sink)
    if verb==1
        progressbar(i,ceil(niter/niter_sink));
    end
    gamma = exp( -(c-r*uq'-up*s')/epsilon );
    b = uq;
    % sinkhorn
    for j=1:niter_sink
        a = p ./ (gamma*b);
        b = q ./ (gamma'*a);
        if nargout>3 || tol>0
            err(end+1) = norm( a.*(gamma*b)-p, 1 );
            if err(end)<tol
                % early stop
                break;
            end
        end
    end
    r = r + epsilon*log(a);
    s = s + epsilon*log(b);
    if ~isempty(err) && err(end)<tol
        % early stop
        break;
    end
end
gamma = exp( -(c-r*uq'-up*s')/epsilon );
err = err(:);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function progressbar(n,N,w)

% progressbar - display a progress bar
%
%    progressbar(n,N,w);
%
% displays the progress of n out of N.
% n should start at 1.
% w is the width of the bar (default w=20).
%
%   Copyright (c) Gabriel Peyr? 2006

if nargin<3
    w = 20;
end

% progress char
cprog = '.';
cprog1 = '*';
% begining char
cbeg = '[';
% ending char
cend = ']';

p = min( floor(n/N*(w+1)), w);

global pprev;
if isempty(pprev)
    pprev = -1;
end

if not(p==pprev)
    ps = repmat(cprog, [1 w]);
    ps(1:p) = cprog1;
    ps = [cbeg ps cend];
    if n>1
        % clear previous string
        fprintf( repmat('\b', [1 length(ps)]) );
    end
    fprintf(ps);
end
pprev = p;
if n==N
    fprintf('\n');
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = getoptions(options, name, v, mandatory)

% getoptions - retrieve options parameter
%
%   v = getoptions(options, 'entry', v0);
% is equivalent to the code:
%   if isfield(options, 'entry')
%       v = options.entry;
%   else
%       v = v0;
%   end
%
%   Copyright (c) 2007 Gabriel Peyre

if nargin<4
    mandatory = 0;
end

if isfield(options, name)
    v = options.(name);
elseif mandatory
    error(['You have to provide options.' name '.']);
end 

end

