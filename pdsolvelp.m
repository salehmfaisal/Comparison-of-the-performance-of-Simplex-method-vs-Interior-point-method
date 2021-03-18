%
% [xstar,optobj,optbasis,nonbasis0,nonbasisu,totaliters,
%  raystar,ystar,wstar,zstar]=
%             pdsolvelp(A,b,c,u,const,maxiters)
%
% Solves an LP of the form
%
%  min const+cx
%      Ax=b
%      0<=x<=u
%
% using the primal-dual interior point method.  This version of the code
% implements scaling and the removal of redundant constraint
% equations.
% 
% If the LP has an optimal solution, the optimal primal solution is 
% returned in xstar, the optimal objective value is in optobj, 
% and the optimal dual solution is returned in ystar, wstar, zstar.
%
% This code returns NaN's if the LP is infeasible, unbounded, or
% the primal-dual method simply fails to converge.
%  
function [optobj,xstar,ystar,wstar,zstar]=...
    pdsolvelp(A,b,c,u,const,tol,maxiters)
%
% Get the problem size.
%
[m,n]=size(A);
%
% Set default parameters if needed.
%
if ((nargin < 6) | isempty(tol))
  tol=1.0e-6;
end
if ((nargin < 6) | isempty(maxiters))
  maxiters=100;
end
%
% First, scale the LP.
%
[As,bs,cs,us,rows,cols]=scalelp(A,b,c,u);
%
% Next, remove any linearly dependent equations from the constraints.  Make
% sure that A is sparse.  
%
[Ar,br,p,rnk]=rreqns(As,bs);
Ar=sparse(Ar);
%
% Output information about the reduced problem.
%
fprintf('rreqns removed %d equations\n',size(As,1)-size(Ar,1));
%
% Get an initial solution.
%
[x0,s0,y0,w0,z0]=initsol(Ar,br,cs,us);
%
% Next, solve the LP.
%
[xs,ss,ysr,ws,zs]=pdpcub(Ar,br,cs,us,const,x0,s0,y0,w0,z0,tol,maxiters);
%
% Undo the effect of rreqns.  We simply assign a dual multiplier of
% 0 to any redundant constraint that was removed from the model.
%
ys=zeros(size(As,1),1);
ys(p)=ysr;
%
% Setup a dummy ray- this isn't returned by pdpcub
%
rays=zeros(length(xs),1);
%
% Unscale the solution.
%
[xstar,raystar,ystar,wstar,zstar]=unscalesoln(rows,cols,xs,rays,ys,ws,zs);
%
% Compute the objective value.
%
optobj=c'*xstar+const;


