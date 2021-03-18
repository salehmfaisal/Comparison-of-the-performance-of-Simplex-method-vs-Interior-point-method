%
% Solve an example netlib LP, with removal of redundant constraints
% and scaling.  These problems were written as minimization
% problems, so we negate the objective function.
%
% Get the problem name and load it in.
%
fname=input('Enter problem file name: ','s');
load(fname);
%
% Start timing.
%
t=tic;
%
% Solve it.  Note that c is given as a column vector in the file,
% so we transpose it.
%
[xstar,optobj,optbasis,nonbasis0,nonbasisu,totaliters,raystar,ystar,wstar,zstar]=solvelp(A,b,-c',u,-const,[],100);
%
% Count the time.
%
simplextime=toc(t);
%
% Output the results.
%
fprintf('Primal Simplex Solution took %.2f seconds\n', simplextime);
fprintf('Optimal Objective = %.8e\n',optobj);
