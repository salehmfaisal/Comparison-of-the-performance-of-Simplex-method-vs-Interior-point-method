function [x,y,z,iter]=pdpc(A,b,c,x0,y0,z0,tol,maxiter)
%
%  [x,y,z,iter]=pdpc(A,b,c,x0,y0,z0,tol,maxiter)
%
%  Implements a primal-dual predictor-corrector method for LP.
%
%  Inputs:
%     A,b,c        The problem, min c'*x, Ax=b, x>=0
%     x0,y0,z0     Initial solution, with x0,z0>0
%     tol          Convergence tolerance (default 1.0e-6.)
%     maxiter      Maximum number of iterations to run (default 100.)
%
%  Outputs:
%     x,y,z        An optimal solution.
%     iter         Iteration count.
%

%
% Get the problem size.
%
[m,n]=size(A);
%
% First, set tol if it wasn't given.
%
if (nargin < 7)
  tol=1.0e-6;
end
%
% Also set maxiter if it wasn't given.
%
if (nargin < 8)
  maxiter=100;
end
%
% Compute the best permutation of A*Zi*X*A'.  Note that this only needs
% to be done once, since the permutation will never change.
% Furthermore, this is indepedent of the actual values of Zi and X, 
% so we can simply use A*I*I*A'=A*A'.
%
M=A*A';
p=symamd(M);
%
% Initialize the loop 
%
x=x0;
y=y0;
z=z0;
deltay=zeros(size(y));
deltax=zeros(size(x));
deltaz=zeros(size(z));
ret=0;
iter=0;
alphap=1;
alphad=1;
while (iter<maxiter)
  if (((norm(A*x-b)/(1+norm(b))) < tol) & ...
      ((norm(A'*y+z-c)/(1+norm(c)) < tol) & ...
      (abs(c'*x-b'*y)/(1+abs(c'*x))) < tol))
    return
  else
%
% Compute mu.
%
  mu=x'*z/n;
%
%  X=diag(x) and Z=diag(z).  While we're here, compute inv(X) and inv(Z).
%  Also compute a vector of all ones.
%
  X=spdiags(x,0,n,n);
  Z=spdiags(z,0,n,n);
  Xi=inv(X);
  Zi=inv(Z);
  e=ones(n,1);
%
% Compute the system matrix.
%
  M=A*Zi*X*A';
%
%  Apply the permutation.
%
  MP=M(p,p);
%
% Factor M
%
  R=chol(MP);
%
% Compute a primal-dual affine step first.
%
  rhs=-M*y+A*X*Zi*c+b-A*x;
  rhs=rhs(p);
%
% Compute the deltay step.
%
  v=R'\rhs;
  deltay(p)=R\v;
%
% Now, compute the deltax and deltaz steps.  
%
  deltax=Zi*X*(A'*(deltay+y)-c);
  deltaz=-Xi*Z*deltax-Z*e;
%
%  Find the maximum alphap
%
  alphap=1;
  for i=1:n,
    if (x(i)+alphap*deltax(i) < 0),
      alphap=-x(i)/deltax(i);
    end;
  end;
%
% Find the maximum alphad
%
  alphad=1;
  for i=1:n,
    if (z(i)+alphad*deltaz(i) < 0),
      alphad=-z(i)/deltaz(i);
    end;
  end;
%
% Compute muplus
%
  muplus=(x+alphap*deltax)'*(z+alphad*deltaz)/n;
%
% Compute sigma.
%
  sigma=min([0.5,muplus/mu])*min((muplus/mu)^2,1);
%
% Now, compute the corrector step.
%
  rhs=-M*y+A*X*Zi*c-sigma*mu*A*Zi*e+b-A*x;
  rhs=rhs(p);
%
% Compute the deltay step.
%
   v=R'\rhs;
   deltay(p)=R\v;
%
% Now, compute the deltax and deltaz steps.  
%
  deltax=Zi*X*(A'*(deltay+y)-c)+sigma*mu*Zi*e;
  deltaz=-Xi*Z*deltax-Z*e+sigma*mu*Xi*e;    
%
%  Find the maximum alphap
%
  alphap=1/0.9995;
  for i=1:n,
    if (x(i)+alphap*deltax(i) < 0),
      alphap=-x(i)/deltax(i);
    end;
  end;
  alphap=0.9995*alphap;
%
% Find the maximum alphad
%
  alphad=1/0.9995;
  for i=1:n,
    if (z(i)+alphad*deltaz(i) < 0),
      alphad=-z(i)/deltaz(i);
    end;
  end;
  alphad=0.9995*alphad;
%
% Take the step.
%
    x=x+alphap*deltax;
    y=y+alphad*deltay;
    z=z+alphad*deltaz;
%
% Output some useful information
%
    iter=iter+1;
    pinfeas=norm(A*x-b)/(1+norm(b));
    dinfeas=norm(A'*y+z-c)/(1+norm(c));
    pobj=c'*x;
    dobj=b'*y;
    fprintf('%03d AP:%.3f AD:%.3f PI:%.2e DI:%.2e PO:%.4e DO:%.4e\n',[iter,alphap,alphad,pinfeas,dinfeas,pobj,dobj]);
  end;
end;
