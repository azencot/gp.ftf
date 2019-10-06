classdef LINFTF < handle
    
    properties (Access='public')
        mesh
        
        tau         % time step
        eps         % epsilon
        Bo          % Bond number
        beta        % Navier slip
        
        evap        % evaporation term
        ce          % evaporation constant
        
        DTS         % dynamic time stepping due to CFL & Orestis conditions
        alpha       %
        ID
        a           % Bo*z - Hv
        b           % Bo*cos(theta) - Tv
        L           % - D*G, D is the divergence, G is the gradient
        
        VA          % vertex areas
        IVA
        ITA         % inverse of triangle areas
        
        st          % simulation time (actual time)
        sti         % computation time per step
        ST          % simulation time (simulated time)
                
        cvk         % current velocity
        R
    end
    
    methods
        
        function ftf = LINFTF( meshname, up, sign )
            
            if nargin < 2; up = 2; end
            if nargin < 3; sign = 1; end;
            ftf.mesh = MESH( meshname, up, sign );
            
            nv = ftf.mesh.nv;
            nf = ftf.mesh.nf;
            va = ftf.mesh.va;
            ta = ftf.mesh.ta;
            
            ftf.ID = speye( nv );
            ftf.VA = spdiags(va,0,nv,nv);
            ftf.IVA = spdiags(1./va,0,nv,nv);
            ftf.ITA = spdiags(1./[ta;ta;ta],0,3*nf,3*nf);
        end
        
        function physical_param( ftf, tau, eps, Bo, beta, ce, DTS, alpha, r )
            
            ftf.tau = tau;
            ftf.eps = eps;
            ftf.Bo = Bo;
            
            ftf.beta = beta;
            ftf.ce = ce;
            
            ftf.evap = 0;
            if ftf.ce ~= 0; ftf.evap = 1; end;
            
            ftf.DTS = 0;
            if nargin > 7 
                ftf.DTS = DTS; 
                ftf.alpha = alpha;
            end;
            
            if nargin < 9; r = 1; end;

            % Given the physical parameters, we can compute a,b and L
            ftf.a = ftf.VA*(ftf.Bo * ftf.mesh.z - ftf.mesh.Hv);
            ftf.b = ftf.VA*(ftf.Bo * ftf.mesh.cos_theta - ftf.mesh.Tv);
            ftf.b = spdiags(ftf.b,0,ftf.mesh.nv,ftf.mesh.nv);
            ftf.L = r*ftf.VA*(- ftf.mesh.D * ftf.mesh.G);
            
            ftf.st = 0;
            ftf.sti = [];
            ftf.ST = [];
        end
        
        function NK = compute_mobility( ftf, u )
            if ftf.beta == 0; u = max(u,1e-3); end;
            
            % standard average
            Nk1 = repmat( ftf.mesh.v2f( ftf.beta + u/3 ), 3, 1 );
            Nk2 = repmat( ftf.mesh.v2f( ftf.eps*(u.^2)/12), 3, 1 );
            Nk2 = ftf.mesh.HS * Nk2;
                        
            NK = ftf.ITA * (Nk1 + Nk2);
            NK = spdiags(NK,0,3*ftf.mesh.nf,3*ftf.mesh.nf);   % 3*nf x 3*nf
        end
        
        function BD = compute_adjoint( ftf, u )
            nv = ftf.mesh.nv;
            
            BD = ftf.mesh.vf2adop( u );
            BD = BD + spdiags( u, 0, nv, nv ) * ftf.mesh.D;
        end
        
        function update_tau( ftf )
            ftf.ST = [ftf.ST ftf.tau];
            if ftf.DTS == 1
                                
                %OV: time stepping from the notes
                cb = diag(ftf.b);
                minb = -min(cb(cb<=0));
                maxva = max(ftf.mesh.va);
                g = norm(ftf.R,'fro')*maxva*(ftf.eps*minb);
               
                tau1 = ftf.alpha/g;
                %OV:end
                
                %OA: CFL condition according to v
                maxv = max( MESH.normv( ftf.cvk ) );
                tau2 = ftf.alpha*ftf.mesh.mel/maxv;
                %OA:end
                
                ftf.tau = min(tau1,tau2);
                
            end
        end
        
        function u = run_sim( ftf, u0, steps )
            
            if ftf.tau < 1e-16 || sum(ftf.VA*u0)<1e-5; u = []; return; end;
            
            outer = tic;

            u = zeros( ftf.mesh.nv, steps );
            u(:,1) = u0;
            
            for k = 1:steps-1
                inner = tic;
                uk = u(:,k);
                
                if ftf.evap == 1
                    uke = exp( - ftf.tau*(uk+ftf.ce).^(-2) ) .* uk;
                else
                    uke = uk;
                end
                
                NK = ftf.compute_mobility( uk );
                BD = ftf.compute_adjoint( uke );
                
                MB = NK*BD';
                ftf.R = BD*MB;

				if ftf.DTS == 1
                    ftf.cvk = MB*(ftf.a + ftf.eps*(ftf.b+ftf.L)*uke);
                end
                
                % linear solve
                u(:,k+1) = ( ftf.ID + ftf.eps*ftf.tau*ftf.R*(ftf.b+ftf.L) ) \ ...
                           ( uke - ftf.tau*ftf.R*ftf.a );
                % dynamic time stepping
                ftf.update_tau();
                
                fprintf('time=%g, step=%d, min(u)=%f, max(u)=%f, int(u)=%f, time step = %e\n', ...
                sum(ftf.ST),k,min(uk),max(uk),sum(ftf.VA*uk),ftf.tau);

                ftf.sti = [ftf.sti toc(inner)];
            
                % break simulation in case tau is too small
                if ftf.tau < 1e-16 || sum(ftf.VA*uk) < 1e-5
                    fprintf('Time step (tau) is smaller then machine EPS\n');
                    u = u(:,1:k);
                    break;
                end
            end
            
            ftf.st = ftf.st + toc(outer);
        end
        
        function [ E ] = compute_energy( ftf, U )
            
            n = size(U,2);
            
            E = zeros(n-1,1);
            for k = 2:n
                
                NK = ftf.compute_mobility( U(:,k-1) );
                BD = ftf.compute_adjoint( U(:,k-1) );
                
                p = ftf.a + ftf.eps*(ftf.b+ftf.L)*U(:,k);
                v = NK * (BD' * p);
                
                E(k-1) = .5*ftf.ST(k-1) * v(:)' * BD' * p + ...
                            ftf.a'*U(:,k) + ...
                         .5*ftf.eps * U(:,k)' * (ftf.b+ftf.L) * U(:,k);
            end
        end
        
        function [ E ] = compute_energy_no_mobility( ftf, U )
            
            n = size(U,2);
            
            E = zeros(n,1);
            for k = 1:n
                E(k) = ftf.a'*U(:,k) + ftf.eps/2*U(:,k)'*ftf.b*U(:,k) + ...
                    ftf.eps/2*U(:,k)'*ftf.L*U(:,k);
            end
        end
        
        function [ P ] = compute_pressure( ftf, U )
            
            n = size(U,2);
            
            P = zeros(ftf.mesh.nv,n);
            for k = 1:n
                P(:,k) = ftf.a + ftf.eps*(ftf.b+ftf.L)*U(:,k);
            end
            P = ftf.IVA*P;
        end
        
        function [ v ] = compute_velocity( ftf, u )
            
            NK = ftf.compute_mobility( u );
            BD = ftf.compute_adjoint( u );

            v = NK*(BD'*(ftf.a + ftf.eps*(ftf.b+ftf.L)*u));
        end
    end
    
end