%
% Demonstrate the importance of sparse Cholesky factorization and
% ordering before factorization in the primal-dual method.
%
%
% First, load in the test problem and determine the size of A.
%
load e226.mat
asize=size(A)
%
% Compute A*A'
%
AAT=A*A';
%
% Compute inv(A*A')
%
t=cputime;
AATINV=inv(AAT);
invtime=cputime-t
%
% Compute the Cholesky factorization of AAT.
%
t=cputime;
R1=chol(AAT);
choltime1=cputime-t
%
% Compute a better ordering of the rows and columns of AAT.
%
t=cputime;
p=symamd(AAT);
AATP=AAT(p,p);
permtime=cputime-t
%
% Compute the Cholesky factorization of AATP.
%
t=cputime;
R2=chol(AATP);
choltime2=cputime-t
%
% Plot spy plots for the various matrices.
%
figure(1);
spy(AAT);
print -deps spyaat.eps
%eps2pdf('spyaat.eps');
figure(2);
spy(AATINV);
print -deps spyaatinv.eps
%eps2pdf('spyaatinv.eps');
figure(3);
spy(R1);
print -deps spyr1.eps
%eps2pdf('spyr1.eps');
figure(4);
spy(AATP);
print -deps spyaatp.eps
%eps2pdf('spyaatp.eps');
figure(5);
spy(R2);
print -deps spyr2.eps
%eps2pdf('spyr2.eps');
%
% Print out the numbers of nonzeros in various matrices.
%
nnzaat=nnz(AAT)
nnzaatinv=nnz(AATINV)
nnzr1=nnz(R1)
nnzR2=nnz(R2)
