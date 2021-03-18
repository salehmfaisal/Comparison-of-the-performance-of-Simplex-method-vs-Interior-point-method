%
% Solve a netlib example problem, both by the primal simplex method
% and the interior point method.
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
[xstar,optobj,optbasis,nonbasis0,nonbasisu,totaliters,ray,ystar,wstar,zstar]=solvelp(A,b,c',u,const,[],100);
%
% Count the time.
%
simplextime=toc(t);
%
% Output information about the simplex solution.
%
fprintf('\n');
fprintf('Min(xstar)=%e\n',min([xstar; 0]));
fprintf('Min(u-xstar)=%e\n',min([u-xstar; 0]));
fprintf('pfeas=%e\n',norm(A*xstar-b)/(1+norm(b)));
fprintf('pobj=%e\n',full(c'*xstar));
%fprintf('dfeas=%e\n',norm(ystar*A-wstar+zstar-c')/(1+norm(c)));
%fprintf('dobj=%e\n',full(ystar*b-wstar*u));
fprintf('Primal Simplex Optimal Objective=%.15e\n',full(optobj));
fprintf('Primal Simplex Solution took %.2f seconds\n', ...
	simplextime);
fprintf('\n');
