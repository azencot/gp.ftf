
paths

meshname = 'wine_glass_tri_r2';
ftf = LINFTF( meshname, 2 );
X = ftf.mesh.vertices; x = X(:,1); z = X(:,2); y = X(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical parameters and initial conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = 1e-8; eps = .0001; Bo = 500; beta = 0; ce = 0; DTS = 1; alpha = 1; r = 100;
ftf.physical_param(tau,eps,Bo,beta,ce,DTS,alpha,r);

vi = find(abs(z - 0.7) < 1e-2);

umin = 0.001; umax = 1; a = .01; 
m = 25; ph = .05;
u = quadratic_surf2( ftf.mesh, vi, umax, umin, a, 200 );
u = u - umin;
ur = (ph/2+ph/5)*sin(m*atan2(y,x)+pi/(m-1))+ ph/2*sin((m-m/2)*atan2(y,x));
u = u .* (1 + ur);
u = u + umin;

%%%%%%%%%%%%%%%%%%
% run simulation %
%%%%%%%%%%%%%%%%%%
fprintf('Computing simulation..\n');

k = 5; steps = 100;

for i = 1:k
    u2 = ftf.run_sim( u(:,end), steps );
    u = [u u2(:,2:end)];
end
figure; MESH_VIS.func(ftf.mesh,u(:,450));

params_str = sprintf('eps_%g_Bo_%g_beta_%g_ce_%g_tau_%g_DTS_%d_alpha_%g_u_%g_r_%g', ...
                      eps, Bo, beta, ce, tau, DTS, alpha, umin, r);
matfilename = sprintf('./experiments/%s_%s.mat', meshname, params_str);
save(matfilename,'u','ftf');