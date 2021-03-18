%
% [x0,s0,y0,w0,z0]=initsol(A,b,c,u)
%
%  Inputs:
%     A,b,c,u,           The problem, min const+c'*x, Ax=b, x>=0, x<=u.
%
%  Outputs:
%     x0,s0,y0,w0,z0     An initial solution.
%
function [x0,s0,y0,w0,z0]=initsol(A,b,c,u)
%
%
[m,n]=size(A);
e=ones(n,1);
z1=max(abs(b'*A*e)/norm(A*e,2),1);
z2=1;

x0=zeros(n,1);
s0=zeros(n,1);
y0=zeros(m,1);
w0=zeros(n,1);
z0=zeros(n,1);

for i=1:n,
  x0(i)=min(z1,u(i)/2);
  s0(i)=u(i)-x0(i);
end;

for i=1:n,
  z0(i)=max(abs(c(i)),z2);
  w0(i)=max(z0(i)-c(i),z2);
end;

for i=1:n,
  if (u(i) > 1.0e30),
    w0(i)=0.0;
    s0(i)=0.0;
  end;
end;
