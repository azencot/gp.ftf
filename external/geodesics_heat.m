function r = geodesics_heat(mesh,vi)

G = mesh.G;     % discrete gradient
D = mesh.D;     % discrete divergence
L = D * G;      % discrete Laplacian (notice the sign)
A = spdiags(mesh.va,0,mesh.nv,mesh.nv);
W = A * L;

t = mesh.mel^2;
u0 = zeros(mesh.nv,1); u0(vi) = 1;

% integrate the heat flow
ut = (A - t*W) \ u0;
gut = reshape(G * ut, mesh.nf, 3);

% normalize and negate
V = - MESH.normalize_vf( gut );

% solve the Poisson equation
b = A * D * V(:);
phi = W \ b;
r = phi-min(phi);     % project out constant

% % exact distance on a sphere 
% r = max(max(abs(X)));
% x = repmat(X(vi,:),mesh.nv,1);
% f = r*acos(dot(x,X,2)/r^2);

% % compute cotan Laplacian
% [W,AA] = cotLaplacian(mesh);
% A = spdiags(AA,0,mesh.nv,mesh.nv);

% % handle boundary
% [ii,ib] = find_bdry(mesh);

% con = zeros(size(ib));

% % Dirichlet BC
% ut1 = solve_poisson_con(A+t*W,u0,ii,ib,con);

% % Neumann BC
% ut2 = solve_poisson_con(A+t*W,u0,ii,[],0,ib,con);