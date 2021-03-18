%
% x=solvebxb(L,U,p,b)
%
% Given an LU factorization L*U=B(p,:), solves Bx=b.
%
function x=solveBxb(L,U,p,b)
w=L\(b(p));
x=U\w;
