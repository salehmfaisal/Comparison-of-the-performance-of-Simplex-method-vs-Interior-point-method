%
% y=solveyBc(L,U,p,c)
%
% Solves the equation y*B=c for y, given the LU factorization L*U=B(p,:)
%
function y=solveyBc(L,U,p,cb)
w=cb/U;
v=w/L;
y(p)=v;
