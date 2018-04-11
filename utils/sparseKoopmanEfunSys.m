function dphi = sparseKoopmanEfunSys(t,y,u,p)
% Delete again?

nvar        = length(y);
polyorder   = p.polyorder;
usesine     = p.usesine;
xi          = p.xi;
phiPool     = buildTheta(y,nvar,polyorder,usesine);
Gamma = buildGamma(x,dx,nvar,polyorder,usesine);
dphi        = (phiPool*xi)';