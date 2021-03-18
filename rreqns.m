%
% [Aout,bout,p,r]=rreqns(A,b);
%
% Input: A system of equations A*x=b.
%
% Output: The system of equations with redundant linearly dependent
% equations removed.  Aout*x=bout. 
%
% The reduced system of equations  A(p,:)*x=b(p) is of full row
% rank.  The function computes Aout=A(p,:) and bout=b(p) for convenience.
%
% The rank is returned in r.
%
% If the system is inconsistent, the routine returns NaN's.
%
% In the context of LP, knowing what equations were kept makes it
% possible to translate y for the reduced problem into y for the
% original problem by yorig(p)=y
%
function [Aout,bout,p,r]=rreqns(A,b);
%
% Get the problem size.
%
[m,n]=size(A);
%
% Test for feasibility of the system of equations.
%
epsilon=1.0e-10;
xtest=A\b;
if (norm(A*xtest-b)/(1+norm(b)) > epsilon)
    warning('Warning A*x=b is inconsistent!');
    Aout=NaN*A;
    bout=NaN*b;
    p=NaN*ones(1,m);
    r=NaN;
    return;
end
%
% Compute the QR factorization.
%
[Q,R,p]=qr(A','vector');
%
% Figure out the rank of R by counting zero rows.
%
r=m;
for i=1:m
  if (norm(R(i,:))==0)
    r=r-1;
  end
end
%
% If r has full rank then just return the original problem.
%
if (r==m)
%
% A has full row rank, so we use the original problem.
%
  Aout=A;
  bout=b;
  p=1:m;
  return
else
%
% A does not have full row rank.
%
  p=p(1:r);
  Aout=A(p,:);
  bout=b(p);
  return
end
