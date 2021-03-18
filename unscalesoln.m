%
% [xstar,raystar]=unscalesoln(r,s,xs,rays,ys,ws,zs)
%
% Unscales the optimal solution to an LP.
%
function [xstar,raystar,ystar,wstar,zstar]=...
    unscalesoln(r,s,xs,rays,ys,ws,zs)
xstar=xs ./ s;
zstar=zs.* s;
wstar=ws.* s;
ystar=ys./ r;
raystar=rays ./ s;