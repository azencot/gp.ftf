function f = quadratic_surf(mesh,vi,a,b,c,d)

n = numel(vi);
% d = 2;
f = 0;
for i = 1:n
    
    r = geodesics_heat( mesh, vi(i) );
    f = f + d*c(i)*max( a/d - r.^2, b/d );
    
end

f = (a-b) * ( (f - min(f)) ./ (max(f)-min(f)) ) + b;