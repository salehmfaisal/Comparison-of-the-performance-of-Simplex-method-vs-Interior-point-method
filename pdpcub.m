%
%  [x,s,y,w,z]=pdpcub(A,b,c,u,const,x0,s0,y0,w0,z0,tol,maxiters)
%
%  Implements a primal-dual predictor correct method.
%
%  Inputs:
%     A,b,c,u,const      The problem, min const+c'*x, Ax=b, x>=0, x<=u.
%     x0,s0,y0,w0,z0     Initial solution, with x0,s0,w0,z0,>0
%                        x+s=u.
%
%  Outputs:
%     x,s,y,w,z          An optimal solution.
%
function [x,s,y,w,z]=pdpcub(A,b,c,u,const,x0,s0,y0,w0,z0,tol, ...
                            maxiters)
%
% Set default values for tol and maxiters.
%
if (nargin < 11)
  tol=1.0e-6;
end
if (nargin < 12)
  maxiters=100;
end
%
% Disable warning about nearly singular matrices.
%
warning('off','MATLAB:SingularMatrix');
warning('off','MATLAB:NearlySingularMatrix');
%
% Initialize some variables.
%
x=x0;
s=s0;
y=y0;
w=w0;
z=z0;
[m,n]=size(A);
ret=0;
iter=0;
alphap=1;
alphad=1;
deltay=zeros(m,1);
deltaw=zeros(n,1);
deltaz=zeros(n,1);
deltax=zeros(n,1);
deltas=zeros(n,1);
%
%  The main loop.
%
while (ret == 0)
%
% First, check for convergence.
%
  if (((norm(A*x-b)/(1+norm(b))) < tol) & ...
      ((norm(A'*y+-w+z-c)/(1+norm(c))) < tol) & ...
      (abs(c'*x-(b'*y-u'*w))/(1+abs(const+c'*x)) < tol))
    ret=1;
    continue;
  else
%
% Check for maxiters
%
    if (iter > maxiters)
      warning('Exceed maxiters');
      x=NaN*x;
      y=NaN*y;
      z=NaN*y;
      s=NaN*s;
      w=NaN*w;
      ret=1;
      continue;
    end
%
% Compute mu.
%
  mu=(x'*z+s'*w)/(2*n);
%
%  Compute the predictor step.
%
  sigma=0; 
%
% Fix up stuff for unbounded variables.
%
  for i=1:n,
    if (u(i) >= 1.0e30),
      s(i)=1.0e30;
      w(i)=0.0;
    end;
  end;
%
%  X=diag(x) and Z=diag(z).  While we're here, compute inv(X) and inv(Z).
%  Also compute a vector of all ones.
%
  X=spdiags(x,0,n,n);
  Z=spdiags(z,0,n,n);
  S=spdiags(s,0,n,n);
  W=spdiags(w,0,n,n);
  
  Xi=inv(X);
  Zi=inv(Z);
  Si=inv(S);
  Wi=inv(W);
  e=ones(n,1);
  rho=sigma*mu*(Si-Xi)*e-(W-Z)*e;
  theta=inv(Xi*Z+Si*W);
%
% Fix up variables with upper bounds of infinity.
%
  for i=1:n,
    if (u(i) >= 1.0e30),
      S(i,i)=1.0e30;
      Si(i,i)=0.0;
      W(i,i)=0.0;
      Wi(i,i)=1.0e30;
      theta(i,i)=x(i)/z(i);
      rho(i)=-sigma*mu/x(i)+z(i);
    end;
  end;

%
% Compute the best permutation of A*Zi*X*A'.  Note that this only needs
% to be done once, since the permutation will never change. 
%
  if (iter == 0),
    p=symamd(A*theta*A');
  end;
%
% Compute the system matrix and rhs.
%
  M=A*theta*A';
  rhs=(b-A*x)+A*theta*(c-A'*y+w-z+rho);
%
%  Apply the permutation.
%
  M=M(p,p);
  rhs=rhs(p);
%
% Compute the deltay step.  If we can't get a Cholesky factor of M,
% then give up.
%
   [R,rnk]=chol(M);
   if (rnk==0)
     foo=R'\rhs;
     deltay(p)=R\foo;
   else
     warning('Cholesky factorization of M failed!');
     ret=1;
     continue;
   end
%
% Now, compute the deltax and deltaz steps.  
%
   deltax=theta*(A'*deltay-rho-c+A'*y-w+z);
   deltas=-deltax;
   deltaz=Xi*(-Z*deltax-X*Z*e+mu*sigma*e);
   deltaw=Si*(-W*deltas-S*W*e+mu*sigma*e);

   for i=1:n,
     if (u(i) >=1.0e30),
       deltaw(i)=0.0;
       deltas(i)=0.0;
     end;
   end;
%
%  Find the maximum alphap
%
    alphap=1;
    for i=1:n,
      if (x(i)+alphap*deltax(i) < 0),
        alphap=-x(i)/deltax(i);
      end;
      if (s(i)+alphap*deltas(i) < 0),
         alphap=-s(i)/deltas(i);
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
      if (w(i)+alphad*deltaw(i) < 0),
        alphad=-w(i)/deltaw(i);
      end;
    end;
%
% Compute mu+
%
    muplus=((x+alphap*deltax)'*(z+alphad*deltaz)+...
            (s+alphap*deltas)'*(w+alphad*deltaw))/(2*n);
%
% Determine sigma.
%
  sigma=min([0.5,muplus/mu])*min((muplus/mu)^2,1);
%
% Now, compute the predictor step.  The left hand side of the
% Newton step equations doesn't change, but the right hand side
% does, so we have to update it.
%
%
% Update the rhs.
%
  rho=sigma*mu*(Si-Xi)*e-(W-Z)*e;
  for i=1:n
    if (u(i)>=1.0e30)
      rho(i)=-sigma*mu/x(i)+z(i);
    end
  end
  rhs=(b-A*x)+A*theta*(c-A'*y+w-z+rho);
  rhs=rhs(p);
%
% Solve for deltay.
%
  foo=R'\rhs;
  deltay(p)=R\foo;
%
% Now, compute the deltax and deltaz steps.  
%
   deltax=theta*(A'*deltay-rho-c+A'*y-w+z);
   deltas=-deltax;
   deltaz=Xi*(-Z*deltax-X*Z*e+mu*sigma*e);
   deltaw=Si*(-W*deltas-S*W*e+mu*sigma*e);

   for i=1:n,
     if (u(i) >=1.0e30),
       deltaw(i)=0.0;
       deltas(i)=0.0;
     end;
   end;
%
%  Find the maximum alphap
%
    alphap=1/0.9995;
    for i=1:n,
      if (x(i)+alphap*deltax(i) < 0),
        alphap=-x(i)/deltax(i);
      end;
      if (s(i)+alphap*deltas(i) < 0),
         alphap=-s(i)/deltas(i);
      end;
    end;
    alphap=alphap*0.9995;
%
% Find the maximum alphad
%
    alphad=1/0.9995;
    for i=1:n,
      if (z(i)+alphad*deltaz(i) < 0),
        alphad=-z(i)/deltaz(i);
      end;
      if (w(i)+alphad*deltaw(i) < 0),
        alphad=-w(i)/deltaw(i);
      end;
    end;
    alphad=alphad*0.9995;

%
% Take the step.
%
    x=x+alphap*deltax;
    s=s+alphap*deltas;
    y=y+alphad*deltay;
    w=w+alphad*deltaw;
    z=z+alphad*deltaz;
%
% Output some useful information
%
    iter=iter+1;
    pinfeas=norm(A*x-b)/(1+norm(b));
    dinfeas=norm(A'*y-w+z-c)/(1+norm(c));
    pobj=const+c'*x;
    dobj=const+b'*y-u'*w;
    fprintf('%03d AP:%.3f AD:%.3f PI:%.2e DI:%.2e PO:%.4e DO:%.4e\n',[iter,alphap,alphad,pinfeas,dinfeas,pobj,dobj]);

  end;
end;
