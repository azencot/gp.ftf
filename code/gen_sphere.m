
paths

meshname = 'sphere_s4';
ftf = LINFTF( meshname, 2 );
X = ftf.mesh.vertices; x = X(:,1); z = X(:,2); y = X(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% physical parameters and initial conditions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tau = 1e-1; eps = 5e-2; Bo = 50; beta = 0; ce = 0; DTS = 1; alpha = 1;
ftf.physical_param(tau,eps,Bo,beta,ce,DTS,alpha);

vi = 13; a = 1; 
umin = .005; umax = .5;
m = 5; ph = 0.5; 
u0str = sprintf('u_%g_m_%g_ph_%g',umin, m, ph);

u = quadratic_surf( ftf.mesh, vi, umax, umin, a, 2 );
u = u - umin;
ur = (ph/2+ph/5)*sin(m*atan2(y,x)+pi/(m-1))+ph/2*sin((m-2)*atan2(y,x));
u = u .* (1 + ur);
u = u + umin;

%%%%%%%%%%%%%%%%%%
% run simulation %
%%%%%%%%%%%%%%%%%%
fprintf('Computing simulation..\n'); 

k = 5; steps = 60;

for i = 1:k
    u2 = ftf.run_sim( u(:,end), steps );
    u = [u u2(:,2:end)];
end
figure; MESH_VIS.func(ftf.mesh,u(:,end));

params_str = sprintf('eps_%g_Bo_%g_beta_%g_ce_%g_tau_%g_DTS_%d_alpha_%g_%s', ... 
                     eps, Bo, beta, ce, tau, DTS, alpha, u0str);
matfilename = sprintf('./experiments/%s_%s.mat', meshname, params_str);
save(matfilename,'u','ftf');