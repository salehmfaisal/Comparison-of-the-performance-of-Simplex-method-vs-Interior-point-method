%
% [xstar,optobj,optbasis,nonbasis0,nonbasisu,totaliters,
%  raystar,ystar,wstar,zstar]=
%             solvelp(A,b,c,u,const,maxiters,printlevel)
%
% Solves an LP of the form
%
%  max const+cx
%      Ax=b
%      0<=x<=u
%
% using the primal simplex method.  This version of the code
% implements scaling and the removal of redundant constraint
% equations.
% 
% If the LP has an optimal solution, the optimal primal solution is 
% returned in xstar, the optimal objective value is in optobj, 
% and the optimal dual solution is returned in ystar, wstar, zstar.
%
% If the LP is unbounded, then optobj=+Inf.  A ray along which the 
% objective value is unbounded is given by xstar+t*raystar.  The 
% dual solution ystar, wstar, zstar is set to NaN.  
%
% If the LP is infeasible, then optobj=NaN.  xtar, raystar, ystar,
% wstar, zstar, optbasis, nonbasis0, and nonbasisu are set to NaN's.
%  
function [xstar,optobj,optbasis,nonbasis0,nonbasisu,totaliters,raystar, ...
          ystar,wstar,zstar]=...
    solvelp(A,b,c,u,const,maxiters,printlevel)
%
% Get the problem size.
%
[m,n]=size(A);
%
% Set default parameters if needed.
%
if ((nargin < 6) | isempty(maxiters))
  maxiters=200*max(m,n/10);
end
if ((nargin < 7) | isempty(printlevel))
  printlevel=0;
end
%
% First, scale the LP.
%
[As,bs,cs,us,r,s]=scalelp(A,b,c,u);
%
% Next, remove any linearly dependent equations from the constraints.  Make
% sure that A is sparse.  
%
[Ar,br,p,rnk]=rreqns(As,bs);
Ar=sparse(Ar);
%
% Output information about the reduced problem.
%
if (printlevel > 0)
  fprintf('rreqns removed %d equations\n',size(As,1)-size(Ar,1));
end
%
% Next, solve the LP.
%
[xs,optobj,optbasis,nonbasis0,nonbasisu,totaliters,rays,ysr,ws,zs]=...
    simplex(Ar,br,cs,us,const,maxiters,printlevel);
%
% Undo the effect of rreqns.  We simply assign a dual multiplier of
% 0 to any redundant constraint that was removed from the model.
%
ys=zeros(1,size(As,1));
ys(p)=ysr;
%
% Unscale the solution.
%
[xstar,raystar,ystar,wstar,zstar]=unscalesoln(r,s,xs,rays,ys,ws,zs);


