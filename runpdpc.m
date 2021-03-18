%
% Solve a netlib example problem, by the primal-dual predictor
% correct method.
%
% This problem requires preprocessing to remove linearly dependent
% columns.  We will also scale the problem.
%
fname=input('Enter problem file name: ','s');
load(fname);
%
% Get problem size
%
[m,n]=size(A);
fprintf('The original problem has %d constraints and %d variables.\n',[m; n]);
%
% Start the clock.
%
t=tic;
%
% Solve it.
%
[optobj,xstar,ystar,wstar,zstar]=pdsolvelp(A,b,c,u,const,[],100);
%
% Count the time.
%
pdtime=toc(t);
%
% Output information about the primal-dual solution.
%
fprintf('\n');
fprintf('Min(xstar)=%e\n',min([xstar; 0]));
fprintf('Min(u-xstar)=%e\n',min([u-xstar; 0]));
fprintf('pfeas=%e\n',norm(A*xstar-b)/(1+norm(b)));
fprintf('pobj=%e\n',full(c'*xstar)+const);
fprintf('dfeas=%e\n',norm(A'*ystar-wstar+zstar-c)/(1+norm(c)));
fprintf('dobj=%e\n',full(ystar'*b-wstar'*u)+const);
fprintf('Primal-Dual Optimal Objective=%.15e\n',full(optobj));
fprintf('Primal-Dual Solution took %.2f seconds\n', ...
	pdtime);
fprintf('\n');
