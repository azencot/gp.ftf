function f = quadratic_surf2(mesh,vi,a,b,c,d)
    
r = geodesics_heat( mesh, vi );
f = d*c*max( a/d - r.^2, b/d );
f = (a-b) * ( (f - min(f)) ./ (max(f)-min(f)) ) + b;