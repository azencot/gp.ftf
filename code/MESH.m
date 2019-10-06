classdef MESH < handle
    
    properties (Access='public')
        
        name
        
        vertices
        triangles
        
        nv              % #vertices
        nf              % #faces
        
        va              % vertex areas
        ta              % face areas
        
        Nf              % normal-per-face
        Nv              % normal-per-vertex
        
        Sf              % shape operator-per-face (vector)
        Sf2             % second shape operator-per-face (vector)
        sf              % shape operator-per-face (tensor)
        sf2             % second shape operator-per-face (tensor)
        Hv              % vertex mean curvature
        Hf              % face mean curvature
        Kv              % vertex Gaussian curvature
        Kf              % face Gaussian curvature
        
        Tv              % vertex total curvature
        HS              % mean curvature + shape operator (tensor)
        z               % up direction
        cos_theta       % Nv \cdot z
        
        G               % tangential gradient operator
        D               % tangential divergence operator
        
        mel             % mean edge length
    end
    
    properties (Access='protected')
        
        E1, E2, E3      % triangle edges

    end
    
    methods
        
        function [ mesh ] = MESH( meshname, up, sign )
            
            if nargin < 1; meshname = 'sphere_s3'; end;
            if nargin < 2; up = 2; end;
            if nargin < 3; sign = 1; end;
            
            mesh.name = meshname;
            
            [mesh.vertices, mesh.triangles] = MESH_READER.readOff([meshname '.off']);
            
            mesh.nv = size(mesh.vertices,1);
            mesh.nf = size(mesh.triangles,1);
            
            mesh.Nf = face_normals( mesh );
            mesh.ta = face_areas( mesh );
            mesh.Nv = vertex_normals( mesh );
            mesh.va = vertex_areas( mesh );
            
            [mesh.E1,mesh.E2,mesh.E3] = face_edges( mesh );
            
            [mesh.Sf,mesh.Sf2] = shape_operator( mesh );
            [mesh.Hv,mesh.Hf] = mean_curvature( mesh );
            [mesh.Kv,mesh.Kf] = gaussian_curvature( mesh );
            [mesh.sf,mesh.sf2] = shapeop_curvature( mesh );
            
            mesh.Tv = total_curvature( mesh );
            mesh.HS = tensor_curvature( mesh );
            [mesh.z,mesh.cos_theta] = gravity( mesh, up, sign );
            
            mesh.G = grad( mesh );
            mesh.D = div( mesh );
        end
        
        function [ ta ] = face_areas( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            ta = mesh.normv( cross( P1, P2 ) ) / 2;
        end
        
        function [ va ] = vertex_areas( mesh )
            va = full( sum( mass_matrix(mesh), 2 ));
        end
        
        function [ M ] = mass_matrix( mesh )
            T = double( mesh.triangles ); 
            
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            Mij = 1/12*[mesh.ta; mesh.ta; mesh.ta];
            Mji = Mij;
            Mii = 1/6*[mesh.ta; mesh.ta; mesh.ta];
            In = [I;J;I];
            Jn = [J;I;I];
            Mn = [Mij;Mji;Mii];
            M = sparse(In,Jn,Mn,mesh.nv,mesh.nv);
        end
                
        function [ Nf ] = face_normals( mesh )
            
            X = mesh.vertices;
            T = mesh.triangles;
            
            P1 = X(T(:,1),:) - X(T(:,2),:);
            P2 = X(T(:,1),:) - X(T(:,3),:);
            
            Nf = cross( P1, P2 );
            Nf = MESH.normalize_vf( Nf );
        end
        
        function [ Nv ] = vertex_normals( mesh )
            TA = spdiags(mesh.ta,0,mesh.nf,mesh.nf);
            
            I = double( repmat(mesh.triangles(:),3,1) );
            J = repmat(1:3,3*mesh.nf,1); J = J(:);
            S = repmat(TA*mesh.Nf,3,1); S = S(:);
            
            Nv = full( sparse(I,J,S,mesh.nv,3) );
            Nv = MESH.normalize_vf( Nv );
        end
        
        function [ E1, E2, E3 ] = face_edges( mesh )
            X = mesh.vertices;
            T = mesh.triangles;
            
            E1 = X(T(:,3),:) - X(T(:,2),:);
            E2 = X(T(:,1),:) - X(T(:,3),:);
            E3 = X(T(:,2),:) - X(T(:,1),:);
            E = [E1; E2; E3];
            
            mesh.mel = mean( MESH.normv( E ) );
        end
        
        function [ Sf, Sf2 ] = shape_operator( mesh )
            if length(mesh.name) >= 6 && strcmp(mesh.name(1:6),'sphere')
                s = 2*eye(3);
                Sf = repmat(s(:),mesh.nf,1);
                Sf2 = zeros(9*mesh.nf,1);
            else
                T = mesh.triangles;
                ITA = repmat(1./mesh.ta,1,3);

                gE = cell(3,1);
                gE{1} = .5 * ITA .* mesh.rotate_vf( mesh.E1 );
                gE{2} = .5 * ITA .* mesh.rotate_vf( mesh.E2 );
                gE{3} = .5 * ITA .* mesh.rotate_vf( mesh.E3 );

                Sf = zeros(9*mesh.nf,1);
                Sf2 = zeros(9*mesh.nf,1);
                for j = 1:mesh.nf
                    s = 0;
                    for i = 1:3
                        s = s + mesh.Nv(T(j,i),:)'*gE{i}(j,:);
                    end
                    % project into the tangent plane
                    p = eye(3) - mesh.Nf(j,:)'*mesh.Nf(j,:);
                    s = p * s;
                    % symmetrize
                    s = (s+s')/2;
                    
                    nx = [0             -mesh.Nf(j,3)    mesh.Nf(j,2); ...
                          mesh.Nf(j,3)  0               -mesh.Nf(j,1); ...
                          -mesh.Nf(j,2) mesh.Nf(j,1)    0];
                    s2 = nx*s*nx;

                    % NOTE: remember this convention!!
                    Sf(9*(j-1)+1:9*j) = - s(:);
                    Sf2(9*(j-1)+1:9*j) = - s2(:);
                end
            end
        end
        
        function [ sf, sf2 ] = shapeop_curvature( mesh )
            I = repmat(1:mesh.nf,3,1)';
            I = [I I+mesh.nf I+2*mesh.nf]';
            I = I(:);

            J = ( 1:mesh.nf )'; 
            J = [J J+mesh.nf J+2*mesh.nf]'; 
            J = repmat(J,3,1);
            J = J(:);
            
            S = mesh.Sf;
            sf = sparse(I,J,S,3*mesh.nf,3*mesh.nf);
            
            S2 = mesh.Sf2;
            sf2 = sparse(I,J,S2,3*mesh.nf,3*mesh.nf);
        end
        
        function [ Hv, Hf ] = mean_curvature( mesh )
            Hf = zeros(mesh.nf,1);
            for j = 1:mesh.nf
                s = mesh.Sf(9*(j-1)+1:9*j);
                s = reshape(s,3,3);
                Hf(j) = trace(s);
            end
            Hf = 2*spdiags([Hf;Hf;Hf],0,3*mesh.nf,3*mesh.nf);
            
            Hv = diag(Hf); 
            Hv = full( Hv(1:mesh.nf) ); 
            Hv = mesh.f2v( Hv );
        end
        
        function [ Kv, Kf ] = gaussian_curvature( mesh )
            Kf = zeros(mesh.nf,1);
            for j = 1:mesh.nf
                s = mesh.Sf(9*(j-1)+1:9*j);
                s = reshape(s,3,3);
                Kf(j) = .5*(trace(s)^2 - trace(s*s));
            end
            Kf = spdiags([Kf;Kf;Kf],0,3*mesh.nf,3*mesh.nf);
            
            Kv = diag(Kf); 
            Kv = full( Kv(1:mesh.nf) ); 
            Kv = mesh.f2v( Kv );
        end
        
        function [ Tv ] = total_curvature( mesh )
            Tv = mesh.Hv.^2 - mesh.Kv;
        end
        
        function [ HS ] = tensor_curvature( mesh )
            HS = 7*mesh.Hf - 3*mesh.sf - 5*mesh.sf2;
        end
        
        function [ z, cos_theta ] = gravity( mesh, up, sign )
            if nargin < 3; sign = 1; end;
            
            z = sign*mesh.vertices(:,up);
            cos_theta = sign*mesh.Nv(:,up);
        end
        
        function [ G ] = grad( mesh )
            % G corresponds to eq. (3.9) in Polygon mesh processing book
            I = repmat(1:mesh.nf,3,1);
            II = [I(:); I(:)+mesh.nf; I(:)+2*mesh.nf];

            J = double( mesh.triangles' );
            JJ = [J(:); J(:); J(:)];

            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );

            S = [ RE1(:) RE2(:) RE3(:) ]';
            SS = S(:);

            G = sparse(II,JJ,SS,3*mesh.nf,mesh.nv);

            TA = .5 * repmat(1 ./ mesh.ta,1,3);
            A = spdiags(TA(:),0,3*mesh.nf,3*mesh.nf);

            G = A*G;
        end
        
        function [ D ] = div( mesh )
            % D corresponds to eq. (3.12) in Polygon mesh processing book
            IVA = spdiags(1./mesh.va,0,mesh.nv,mesh.nv);
            TAC = spdiags(repmat(mesh.ta,3,1),0,3*mesh.nf,3*mesh.nf);

            D = - IVA * mesh.G' * TAC;
        end
        
        function [ op ] = vf2op( mesh, vf )
            vf = reshape( vf, mesh.nf, 3 );
            
            RE1 = mesh.rotate_vf( mesh.E1 );
            RE2 = mesh.rotate_vf( mesh.E2 );
            RE3 = mesh.rotate_vf( mesh.E3 );
            
            T = mesh.triangles;
            I = [T(:,1);T(:,2);T(:,3)];
            J = [T(:,2);T(:,3);T(:,1)];
            Sij = 1/6*[dot(RE2,vf,2); dot(RE3,vf,2); dot(RE1,vf,2)];
            Sji = 1/6*[dot(RE1,vf,2); dot(RE2,vf,2); dot(RE3,vf,2)];
            In = [I;J;I;J];
            Jn = [J;I;I;J];
            Sn = [Sij;Sji;-Sij;-Sji];
            W = sparse(In,Jn,Sn,mesh.nv,mesh.nv);
            IVA = spdiags(1./mesh.va,0,mesh.nv,mesh.nv);
            
            op = IVA*W;
        end
        
        function [ adop ] = vf2adop( mesh, f )
            % the following holds:
            % D_v*u = adop_u*v, adop \in nv x 3*nf
            
            g = mesh.G * f;

            I = repmat(1:mesh.nf,1,3)'; I = I(:);
            J = ( 1:3*mesh.nf )';
            S = g(:);
            BD = sparse(I,J,S,mesh.nf,3*mesh.nf);

            I = double( mesh.triangles(:) );
            J = 1:mesh.nf; J = repmat(J(:),3,1);
            S = repmat(mesh.ta,3,1);
            IVA = spdiags(1./(3*mesh.va),0,mesh.nv,mesh.nv);
            P = IVA*sparse(I,J,S,mesh.nv,mesh.nf);
            
            adop = P * BD;
        end
        
        function [ ff ] = v2f( mesh, fv )
            ff = mean( fv(mesh.triangles), 2 );
        end
    
        function [ fv ] = f2v( mesh, ff )
            F = repmat(ff.*mesh.ta/3,1,3);

            T = double( mesh.triangles );
            I = [T(:,1); T(:,2); T(:,3)];
            J = ones(size(I));
            S = [F(:,1); F(:,2); F(:,3)];

            fv = sparse(I,J,S,mesh.nv,1);
            fv = spdiags(1./mesh.va,0,mesh.nv,mesh.nv)*fv;
            fv = full(fv);
        end
        
        function [ rvf ] = rotate_vf( mesh, vf )
            rvf = cross( mesh.Nf, vf );
        end
        
    end
    
    methods (Static)
        
        function [ nv ] = normv( v )
            nv = sqrt(sum(v.^2,2));
        end
        
        function [ nnv ] = normalize_vf( v )
            nnv = v ./ repmat(MESH.normv(v),1,3);
        end

    end
end