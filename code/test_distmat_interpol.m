%%
% Test for the computation of interpolation of several distances matrices.
% Tests for both 3-D meshes and 2-D point clouds as inputs.

addpath('toolbox/');
addpath('mexEMD/');
addpath('mesh2d/');
addpath('toolbox_geometry/');
addpath('data/meshes/');
addpath('data/meshes-nonrigid/');
addpath('data/shapes/');
addpath('data/shapes-nonrigid/');

%%
% Main parameters.

if not(exist('names'))
    names = {'centaur1' 'horse1'};
end
if not(exist('dist_mode'))
    dist_mode = 'geodesic';
    dist_mode = 'euclidean';
end
if not(exist('init')) % which cloud is used for initialization
    init = 1;
end
% loss being used 
if not(exist('gw_loss'))
    gw_loss = 'kl';
    gw_loss = 'l2';
    sigma = .8;
end
if strcmp(gw_loss(1:2),'kl') 
    if not(isempty(gw_loss(3:end)))
        sigma = str2num(gw_loss(3:end))/10;
    end
    gw_loss = 'kl';
end
% embedding method
mdsmethod = 'classical';
mdsmethod = 'smacof';

%% 
% rep creation

Nshape = length(names);
name = '';
for i=1:Nshape
	name = [name names{i}]; 
    if not(i==Nshape)
        name = [name '-'];
    end
end

lossstr = gw_loss;
if strcmp(gw_loss,'kl')==1
    lossstr = [gw_loss num2str(round(sigma*10))];
end

repbase = ['results/distmat-interpolation/2d-' dist_mode '/'];

rep = [repbase name '/'];
repsubfig = [rep lossstr '/'];
if not(exist(repsubfig))
    mkdir(repsubfig);
end

%%
% Other parameters.

% size of the input clouds
N0 = 400; 
N0 = 1000; 
% parameter for geodesic distances computations.
t_Varadhan = 1e-2;  % time step for Varadhan formula

options.gw_loss = gw_loss;
% switch from matrix distance to kernel
switch gw_loss
    case 'l2'
        Kembed = @(d)d;
        iKembed = @(K)K;
    case 'kl'
        if not(exist('sigma'))
            sigma = 1;
        end
        Kembed = @(d)exp(-(d/sigma).^2);
        iKembed = @(K)sqrt(max(-log(K),0))*sigma;
    otherwise 
        error('Unknorn loss.');
end

%%
% Helpers.

dotp = @(x,y)sum(x(:).*y(:));
mynorm = @(x)norm(x(:));
normalize = @(x)x/sum(x(:));

%%
% Load shapes.

X = {}; F = {}; D = {}; Mu = {}; N = [];
for i=1:Nshape
    opt.upsampling = 10;
    opt.smooth_iter = 0; % control smoothness of the boundary.
    switch dist_mode
        case 'euclidean'
            opt.meshing_mode = 'farthestpoint'; % crude meshing
        case 'geodesic'
            opt.meshing_mode = 'mesh2d'; % adaptive meshing
            opt.Nbound = 500; % #boundary point
    end    
    [V{i},F{i}] = load_shape_triangulation(names{i}, N0*5, opt);    
    switch dist_mode
        case 'euclidean'
            % use Euclidean distance matrices
            Dfull{i} = distmat(V{i}');
        case 'geodesic'
            % use geodesic distance matrix
            Op = load_mesh_operators(V{i}',F{i});
            U = inv( full( diag(Op.AreaV) + t_Varadhan * Op.Delta ));
            U = -t_Varadhan*log(U);
            Dfull{i} = sqrt(max(U,0));
    end    
    % perform subsampling
    Isub{i} = 1:size(Dfull{i},1);
    if size(Dfull{i},1)>N0
        % using euclidean distances for subsampling seems better
        Isub{i} = perform_fartherpoint_subsampling_euclidean(V{i},N0,1);
    end
    X{i} = V{i}(Isub{i},:);
    D{i} = Dfull{i}(Isub{i},Isub{i});
    D{i} = D{i}/median(D{i}(:));
    % corresponding kernel
    K{i} = Kembed(D{i});
    %
    N(i) = size(X{i},1);
    Mu{i} = normalize(ones(N(i),1));
end

%%
% Display as point clouds.

for i=1:Nshape
    clf;
    plot_pointcloud(X{i}, X{i});
    axis tight; axis equal; axis off;
    saveas(gcf, [rep name '-input-' num2str(i) '.eps'], 'epsc');
end

%%
% Display distance function.

clf;
for i=1:Nshape
    subplot(1,Nshape,i); hold on;
    [~,k] = min( V{i}(:,2) ); [~,k] = max( Dfull{i}(k,:) ); % find cool seed point
    opt.face_vertex_color = rescale( Dfull{i}(k,:)' );
    opt.view_param = [-50,5];
    plot_mesh(V{i},F{i}, opt);
    if size(V{i},2)==2 % for 2D plots
        view(2); axis ij;
        plot( V{i}(k,1),V{i}(k,2), 'r.', 'MarkerSize', 30 );
    else
        plot3( V{i}(k,1),V{i}(k,2),V{i}(k,3), 'r.', 'MarkerSize', 30 );
    end
    lighting none; colormap parula(256);
end
if 0
    saveas(gcf, [rep name '-' dist_mode '-input-dist.png'], 'png');
end

%%
% Test for one GW computation.

% N=400
if strcmp(gw_loss, 'l2')
    epsilon = .1/14; % 2-D shapes L2, N=400
    epsilon = .1/20;
    epsilon = .1/15;
else
    epsilon = .1/40; % 2-D shapes KL
    epsilon = .1/10; % sigma=.8, ok
    epsilon = .1/5;
end

% epsilon = 0; % exact-EMD
options.niter = 40; 
options.niter_sinkhorn = 100;
options.tol_sinkhorn = 1e-6; % tolerance on marginals for Sinkhorn
options.tol_gw = 1e-4; % early stopping for GW
options.gamma_init = [];
options.verb = 1;
options.log_domain = 0; % use stabilization or not
[gamma,a,b,GW] = perform_gw_sinkhorn(Kembed(D{1}),Kembed(D{2}),Mu{1},Mu{2},epsilon*0, options);

clf;
plot(GW); axis tight;

%%
% Barycenter computation. 

% time for interpolation
tlist = .5;
tlist = 0;
tlist = linspace(0,1,8);
% option for barycenter
mu = normalize( ones(N(init),1) );
d0 = D{init}(:);
DbaryList = {};
options.verb = 0;
options.bary_method = 'alterating';
options.niter_alternating = 10;
options.tol_alternating = 1e-4;
% options for display
mds_init = 'interp'; % use the linear interpolant as initialization
mds_init = 'input'; % use the first cloud as initialization
mycol = @(t)[1-t t 0];

for k=1:length(tlist)
    t = tlist(k);
    lambda = [1-t t];
    fprintf('Barycenter %d/%d ', k, length(tlist));
    [Kbary,R] = compute_gw_barycenters(K,K{init},lambda,epsilon,options);
    Dbary = iKembed(Kbary);
    DbaryList{k} = Dbary;
    % run MDS
    opts.niter = 80; % for smacof
    opts.verb = 1;
   	opts.X0 = X{init};
    fprintf('Rendering  %d/%d ', k, length(tlist));
    [X1,s] = perform_mds(DbaryList{k},size(X{1},2), mdsmethod, opts);
    % display
	clf;
    plot_pointcloud(X1,X{init});    
	saveas(gcf, [repsubfig num2str(k) '.eps'], 'epsc');
end
