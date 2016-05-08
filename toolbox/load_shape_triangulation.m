function [X,F] = load_shape_triangulation(name, N, options)

% load_shape_triangulation - load a shape from a binary image and triangulate it.
%
%   X is of size (N,2)
%   F is of size(3,P)
%
%   set options.smooth_iter to perform iterate mean curvature motion
%   smoothing.
%
%   Copyright (c) 2015 Gabriel Peyre

% HELPERS for curve management
real2cpx = @(c)c(1,:)' + 1i*c(2,:)';
cpx2real = @(gamma)transpose( [real(gamma(:)),imag(gamma(:))] );
curvabs = @(gamma)[0;cumsum( 1e-5 + abs(gamma(1:end-1)-gamma(2:end)) )];
resample1 = @(gamma,d,P)interp1(d/d(end),gamma,(0:P-1)'/P, 'linear');
resample = @(gamma,P)resample1( [gamma;gamma(1)], curvabs( [gamma;gamma(1)] ),P );
% filtering and convolution
cconv = @(a,b)(ifft( fft(a).*fft(b) ));
normalize = @(h)h/sum(h(:));
tlist = @(N)([0:ceil(N/2)-1, -floor(N/2):-1]')/N;
gaussian = @(mu,N)normalize( exp( -(tlist(N).^2)/(2*mu^2) ) );
gaussfilt = @(f,mu)cconv(f(:),gaussian(mu,length(f)));

options.null = 0;
r = getoptions(options, 'upsampling', 5); % upsampling factor
P = getoptions(options, 'img_resol', 512*2);
smooth_iter = getoptions(options, 'smooth_iter', 20);
meshing_mode = getoptions(options, 'meshing_mode', 'farthestpoint');

%%
% Load the image.

f = load_image(name, P);
f = sum(f,3); % convert to grayscale

thresh = @(f)double(f>mean(f(:)));
f = thresh(f);
if f(1)==1
    f = 1-f;
end

switch meshing_mode
    case 'mesh2d'
        
        %%% ADAPTIVE MESHING %%%
        % #boundary points
        Nbound = getoptions(options, 'Nbound', 100);
        % extract boundary
        t = (0:size(f,1)-1)'/size(f,1);
        c = extract_levelset(t, f, 0);
        c = resample(real2cpx(c),length(c));
        % smooth a bit
        mu = .01/5;
        c = gaussfilt(c,mu);
        % resample
        c = cpx2real(resample(c,Nbound));
        % do the meshing
        hdata.hmax = 10; % largest possible inside
        hdata.hmax = norm(c(:,1)-c(:,2))*1.1;
        options.output = false;
        [X,F] = mesh2d(c',[],hdata,options);
        F = F';
        
    case 'farthestpoint'
        
        %%
        % Smooth it a bit.
        
        n = 41; % width of the convolution kernel
        s = 2;
        t = -n:n;
        g = exp(-t.^2 / (2*s^2) ); g = g/sum(g(:));
        xi = @(x)conv2(conv2(x, g, 'same')', g, 'same')';
        for k=1:smooth_iter
            f = thresh(xi(f));
        end
        
        %%
        % Sample.
        
        I = find(f==1); I = I(randperm(length(I))); I = I(1:r*N);
        x = linspace(0,1,P);
        [x,y] = meshgrid(x,x);
        X = [x(I), y(I)];
        % farthestpoint sampling
        I = perform_fartherpoint_subsampling_euclidean(X,N,1);
        X = X(I,:);
        
        %%
        % Build a delaunay triangulation of the inside.
        
        F = delaunay(X(:,1),X(:,2));
        a = ( X(F(:,1),:) + X(F(:,2),:) + X(F(:,3),:) )/3;
        v = interp2(x,y,f,a(:,1),a(:,2));
        F = F(v>.5,:);
        F = F';
        
    otherwise
        
        error('Unknown meshing mode.');
        
end


% normalize
X = X';
X = X - repmat(mean(X,2), [1 size(X,2)]);
X = X ./ std(X(:)); X = X';



end