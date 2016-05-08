function Op = load_mesh_operators(X,F)

% load_mesh_operators - grad/div/laplacian on mesh
%
%   Op = load_mesh_operators(X,F);
%
%   Warning: Op.Grad is actually -Grad. 
%       op.Delta is the symmetric un-normalized Laplacian. 
%   To obtain a "faithful" non symmetric Laplacian, you need to use
%       spdiags(1./op.AreaV,0,N,N) * op.Delta
%
%   Copyright (c) 2016 Gabriel Peyre


n = size(X,2);
m = size(F,2);

if size(X,1)==2
    X = cat(1,X,zeros(1,n));
end

%%
% Callback to get the coordinates of all the vertex of index \(i=1,2,3\) in
% all faces.

XF = @(i)X(:,F(i,:));

%%
% Compute un-normalized normal through the formula 
% \(e_1 \wedge e_2 \) where \(e_i\) are the edges.

Na = cross( XF(2)-XF(1), XF(3)-XF(1) );

%%
% Compute the area of each face as half the norm of the cross product.

amplitude = @(X)sqrt( sum( X.^2 ) );
AreaF = amplitude(Na)/2;

%%
% Compute the set of unit-norm normals to each face.

normalize = @(X)X ./ repmat(amplitude(X), [3 1]);
Normals = normalize(Na);

%%
% Populate the sparse entries of the matrices for the 
% operator implementing \( \sum_{i \in f} u_i (N_f \wedge e_i) \).

I = []; J = []; V = []; % indexes to build the sparse matrices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % vector N_f^e_i
    wi = cross(XF(t)-XF(s),Normals);    
    % update the index listing
    I = [I, 1:m];
    J = [J, F(i,:)]; 
    V = [V, wi];
end

%% 
% Sparse matrix with entries \(1/(2A_f)\).

dA = spdiags(1./(2*AreaF(:)),0,m,m); 

%%
% Compute gradient.

GradMat = {};
for k=1:3
    GradMat{k} = dA*sparse(I,J,V(k,:),m,n);
end

%%
% \(\nabla\) gradient operator.

Grad = @(u)[GradMat{1}*u, GradMat{2}*u, GradMat{3}*u]';

%%
% Compute divergence matrices as transposed of grad for the face area inner product.

dAf = spdiags(2*AreaF(:),0,m,m); 
DivMat = {GradMat{1}'*dAf, GradMat{2}'*dAf, GradMat{3}'*dAf};

%%
% Div operator. 

Div = @(q)DivMat{1}*q(1,:)' + DivMat{2}*q(2,:)' + DivMat{3}*q(3,:)';

%%
% Laplacian operator as the composition of grad and div. 

Delta = DivMat{1}*GradMat{1} + DivMat{2}*GradMat{2} + DivMat{3}*GradMat{3};

%%
% Cotan of an angle between two vectors.

cota = @(a,b)cot( acos( dot(normalize(a),normalize(b)) ) );

%%
% Compute cotan weights Laplacian. 

I = []; J = []; V = []; % indexes to build the sparse matrices
Ia = []; AreaV = []; % area of vertices
for i=1:3
    % opposite edge e_i indexes
    s = mod(i,3)+1;
    t = mod(i+1,3)+1;
    % adjacent edge
    ctheta = cota(XF(s)-XF(i), XF(t)-XF(i));
    % ctheta = max(ctheta, 1e-2); % avoid degeneracy
    % update the index listing
    I = [I, F(s,:), F(t,:)];
    J = [J, F(t,:), F(s,:)]; 
    V = [V, ctheta, ctheta];
    % update the diagonal with area of face around vertices
    Ia = [Ia, F(i,:)];
    AreaV = [AreaV, AreaF]; 
end
% Aread diagonal matrix
Ac = sparse(Ia,Ia,AreaV,n,n);
% Cotan weights
Wc = sparse(I,J,V,n,n);
% Laplacian with cotan weights.
DeltaCot = spdiags(full(sum(Wc))', 0, n,n) - Wc;

% another implementaiton, that compute correctly area weights.
[DeltaCot,A] = cotLaplacian(X',F');

Op.Delta = Delta;
Op.Grad = Grad;
Op.Div = Div;
Op.DivMat = DivMat;
Op.GradMat = GradMat;
Op.Normals = Normals;
Op.AreaF = AreaF(:); % face area
Op.AreaV = full(A); % vertex area
Op.DeltaCot = DeltaCot;

end

function [W,A] = cotLaplacian(X, T)

% Compute the cotangent weight Laplacian.
% W is the symmetric cot Laplacian, and A are the area weights 

nv = size(X,1);
nf = size(T,1);

if size(X,2)==2
    X = cat(2,X,zeros(nv,1));
end

% Find orig edge lengths and angles
L1 = normv(X(T(:,2),:)-X(T(:,3),:));
L2 = normv(X(T(:,1),:)-X(T(:,3),:));
L3 = normv(X(T(:,1),:)-X(T(:,2),:));
EL = [L1,L2,L3];
A1 = (L2.^2 + L3.^2 - L1.^2) ./ (2.*L2.*L3);
A2 = (L1.^2 + L3.^2 - L2.^2) ./ (2.*L1.*L3);
A3 = (L1.^2 + L2.^2 - L3.^2) ./ (2.*L1.*L2);
A = [A1,A2,A3];
A = acos(A);

% The Cot Laplacian 
I = [T(:,1);T(:,2);T(:,3)];
J = [T(:,2);T(:,3);T(:,1)];
S = 0.5*cot([A(:,3);A(:,1);A(:,2)]);
In = [I;J;I;J];
Jn = [J;I;I;J];
Sn = [-S;-S;S;S];

% Compute the areas. Use mixed weights Voronoi areas
cA = 0.5*cot(A);
vp1 = [2,3,1]; vp2 = [3,1,2];
At = 1/4 * (EL(:,vp1).^2 .* cA(:,vp1) + EL(:,vp2).^2 .* cA(:,vp2));

% Triangle areas
N = cross(X(T(:,1),:)-X(T(:,2),:), X(T(:,1),:) - X(T(:,3),:));
Ar = normv(N);

% Use barycentric area when cot is negative
locs = find(cA(:,1) < 0);
At(locs,1) = Ar(locs)/4; At(locs,2) = Ar(locs)/8; At(locs,3) = Ar(locs)/8;
locs = find(cA(:,2) < 0);
At(locs,1) = Ar(locs)/8; At(locs,2) = Ar(locs)/4; At(locs,3) = Ar(locs)/8;
locs = find(cA(:,3) < 0);
At(locs,1) = Ar(locs)/8; At(locs,2) = Ar(locs)/8; At(locs,3) = Ar(locs)/4;

% Vertex areas = sum triangles nearby
I = [T(:,1);T(:,2);T(:,3)];
J = ones(size(I));
S = [At(:,1);At(:,2);At(:,3)];

if strcmp(version,'8.3.0.532 (R2014a)') || strcmp(version,'8.4.0.150421 (R2014b)')
    W = sparse(double(In),double(Jn),Sn,nv,nv);
    A = sparse(double(I),double(J),S,nv,1);
else
    W = sparse(In,Jn,Sn,nv,nv);
    A = sparse(I,J,S,nv,1);
end

end

function nn = normv(V)

nn = sqrt(sum(V.^2,2));

end