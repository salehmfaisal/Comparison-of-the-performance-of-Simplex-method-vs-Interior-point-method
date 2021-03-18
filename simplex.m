%
% [x,optobj,optbasis,nonbasis0,nonbasisu,totaliters,ray,y,w,z]=
%       simplex(A,b,c,u,const,maxiters,printlevel)
%
% Solves an LP of the form 
%
% min c*x+const
%     Ax=b
%      x >= 0
%      x <= u          (upper bounds may be +Inf)  
%
% by the two-phase simplex method.
%
% Inputs:
%   A,b,c,u,const     Problem data.
%   maxiters          Maximum number of iterations, defaults to
%                     10*m for each phase. 
%   printlevel        if 0, output nothing at all.
%                     if >0, output information every printlevel 
%                     iterations and at the start and end of each phase.
%                     e.g. "400: w=27.2" if w is 27.2 at the
%                     completion of the 400th iteration of phase I.
%                     Defaults to 0.
%
% Outputs:
%  1. If the LP has an optimal BFS,
%     x              optimal solution.
%     z              optimal objective value.
%     optbasis       row vector listing basic variables.
%     totaliter      total phase I and phase II iterations.
%     ray            zeros(n,1)
%
%  2. If the LP is infeasible,
%     x              NaN*ones(n,1)
%     z              NaN
%     optbasis       []
%     totaliters     Phase I iterations.
%     ray            NaN*ones(n,1)
%
%  3. If the LP is unbounded,
%     x              Last BFS.
%     z              -Inf
%     optbasis       Last basis
%     totaliters     total phase I and phase II iterations.
%     ray            optimal ray.  c*(x+t*ray) -> -infinity as
%                    t->infinity
%
%  4. If max iterations exceeded in either phase, exit with
%     error()
%
function [x,optobj,optbasis,nonbasis0,nonbasisu,totaliters,ray,y,w,z]=simplex(A,b,c,u,const,maxiters,printlevel)
%
% Zero tolerances.
%
epsilon1=1.0e-6;
epsilon2=1.0e-6;
epsilon3=1.0e-6;
%
% Get basic problem size data.
%
[m,n]=size(A);
%
%%PhaseI Lp
%% 
Ap = auxi(A,b); %updates A with auxilliary variables for phase I 
maxiters = 10*m;
basis = n+1:n+m;
nonbasis0 = 1:n;
nonbasisu = [];
B = Ap(:,basis);
AN0 = Ap(:,nonbasis0);
ANU = Ap(:,nonbasisu);
c1 = zeros(1,n);
c2 = ones(1,m);
cp = [c1 c2];
u1 = Inf(m,1);
up = vertcat(u,u1);
ub = up(basis);
cb = cp(basis);
cn0 = cp(nonbasis0);
cnu = cp(nonbasisu);
[L,U,p] = lu(B,'vector');
y = solveyBc(L,U,p,cb);
xu = u(nonbasisu);
if isempty(ANU) == 1
    xb = solveBxb(L,U,p,b);
else
    nb = b - ANU*xu;
    xb = solveBxb(L,U,p,nb);
end
rn0 = cn0 - y*AN0;
rnu = [];
if isempty(nonbasisu) == 1
else
    rnu = cnu - y * ANU;
end
mr0 = min(rn0);
mru = max(rnu);
j = 0;
fprintf('Starting Phase I\n');
 xp = zeros(n+m,1);
 xp(basis)= xb;
 xu = u(nonbasisu);
xp(nonbasisu)= xu;
 w = cp*xp;

   fprintf('Iteration =%d, w= %d\n', j, w);

while mr0 < -epsilon1 || (~isempty(nonbasisu) && mru > epsilon1)
if max(xb - ub) < epsilon3 && min(xb) > -1*epsilon3
[~, opt0] = min(rn0); %most negative coefficient
[~, optu] = max(rnu);
if abs(mr0) > abs(mru)
    enteringvar = nonbasis0(opt0);
     a = Ap(:,enteringvar); 
    d = solveBxb(L,U,p,a);
elseif abs(mru) > abs(mr0)
    enteringvar = nonbasisu(optu);
    a = Ap(:,enteringvar); 
    d = solveBxb(L,U,p,a);
    d = -1*d;
else
    enteringvar = nonbasis0(opt0);
     a = Ap(:,enteringvar); 
    d = solveBxb(L,U,p,a);
end
limitentering = up(enteringvar);
[tlimit,leavingvar,leavingbound]=ratiotestub(basis,xb,ub,d,limitentering,epsilon2,epsilon3);
 basis = setdiff(basis, leavingvar);
nonbasis0 = setdiff(nonbasis0, enteringvar);
nonbasisu = setdiff(nonbasisu, enteringvar);
 if ismember(enteringvar, nonbasis0)==1
nonbasis0 = setdiff(nonbasis0, enteringvar);
else 
nonbasisu = setdiff(nonbasisu, enteringvar);
end
if leavingvar == 0
    if leavingbound == 0
       nonbasisu = [nonbasisu enteringvar];
    else
      nonbasis0 = [nonbasis0 enteringvar];
     end
        
else
basis = setdiff(basis, leavingvar);
basis = [basis enteringvar];

 if leavingbound == 0
      nonbasis0 = [nonbasis0 leavingvar];
 else
     nonbasisu = [nonbasisu leavingvar];
 end
 
end
 B = Ap(:,basis);
AN0 = Ap(:,nonbasis0);
ANU = Ap(:,nonbasisu);
cb = cp(basis);
cn0 = cp(nonbasis0);
cnu = cp(nonbasisu);
ub = up(basis);
[L,U,p] = lu(B,'vector');
y = solveyBc(L,U,p,cb);
rnu = cnu - y*ANU;
rn0 = cn0 - y*AN0;
xu = u(nonbasisu);
if isempty(ANU) == 1
    xb = solveBxb(L,U,p,b);
else
    nb = b - ANU*xu;
    xb = solveBxb(L,U,p,nb);
end 
 mr0 = min(rn0);
 mru = max(rnu);
 xp = zeros(n+m,1);
xp(basis)= xb;
xu = u(nonbasisu);
xp(nonbasisu)= xu;
 w = cp*xp;
j = j+1;
if mod(j,printlevel) == 0
   level=sprintf('Iteration =%d, w= %d',...
      j, w);
   disp(level);
end
end
end
w = cp*xp;
if w ~= 0 % Infeasible LP
    optbasis = [];
    x = NaN*ones(n,1);
    z = NaN;
    totaliters = j;
    ray = NaN*ones(n,1);
    optobj = z; %change this
    fprintf('w=%d,Infeasible solution for phase1\n',w);
end 

fprintf('Phase I finished in %d iterations\n',j);

% 
[auv, ii] = max(basis); %Removing auxiliary variable from the from the basis of an optimal solution
while auv > n & w == 0
    for i = 1:n
        a10= ismember(i,nonbasis0); 
        alu = ismember(i,nonbasisu);
        if a10 == 1 || alu == 1
           enteringvar = i;
           a = Ap(:,enteringvar);
           d = B\a;
        end
       if abs(d(ii)) > epsilon1
            basis = setdiff(basis, auv);
            nonbasis0 = setdiff(nonbasis0, enteringvar);
            nonbasisu = setdiff(nonbasisu, enteringvar);
            basis = [basis enteringvar];
           nonbasis0 = [nonbasis0 auv];
           break
       end
    end
 B = Ap(:,basis);
[auv, ii] = max(basis);

end



%% Phase II
% Set default parameters if needed.
%
if ((nargin < 6) | isempty(maxiters))
  maxiters=10*m;
end
if ((nargin < 7) | isempty(printlevel))
  printlevel=0;
end

if  w == 0   % Optimal solution for phase 1
    fprintf('Starting Phase II\n');
    nonbasis0 = setdiff(nonbasis0,n+1:n+m);
    nonbasisu = setdiff(nonbasisu,n+1:n+m);
    B = A(:,basis);
    ANU = A(:,nonbasisu);
    AN0 = A(:,nonbasis0);
    cb = c(basis);
    cn0 = c(nonbasis0);
    cnu = c(nonbasisu);
    ub = u(basis);
   [L,U,p] = lu(B,'vector');
   y = solveyBc(L,U,p,cb);
       xu = u(nonbasisu);
    if isempty(ANU) == 1
        xb = solveBxb(L,U,p,b);
    else
        nb = b - ANU*xu;
        xb = solveBxb(L,U,p,nb);
    end
   rn0 = cn0 - y*AN0;
   if isempty(nonbasisu) == 0
    rnu = cnu - y*ANU;
   end
   mr0 = min(rn0);
   mru = max(rnu);
   k = 0;
   x = zeros(n,1);
   x(basis)= xb;
   xu = u(nonbasisu);
   x(nonbasisu)= xu;
   z = c*x+const;
   fprintf('Iteration =%d, z= %d\n',k, z);

while mr0 < -epsilon1 || (~isempty(mru) && mru > epsilon1)
if max(xb - ub) < epsilon3 && min(xb) > -epsilon3
[~, opt0] = min(rn0); %most negative coefficient
[~, optu] = max(rnu);
if abs(mru) > abs(mr0)
    enteringvar = nonbasisu(optu);
     a = A(:,enteringvar); 
    d = solveBxb(L,U,p,a);
    d = -1*d;
else 
    enteringvar = nonbasis0(opt0);
    a = A(:,enteringvar); 
    d = solveBxb(L,U,p,a);
end

limitentering = u(enteringvar);
[tlimit,leavingvar,leavingbound]=ratiotestub(basis,xb,ub,d,limitentering,epsilon2,epsilon3);
if tlimit == +Inf;
    x = zeros(n,1);
    x(basis)=xb;
    xu = u(nonbasisu);
    x(nonbasisu)= xu;
    z = -Inf;
    t = enteringvar;
    basis = sort(basis);
    optbasis = basis;
    ray = zeros(n,1);
    ray1 = rref(A);
    ray(t) = 1;
    ray(basis) = -ray1(:,t);
     break;
end

if leavingvar == 0
    if leavingbound == 1
    if ismember(enteringvar, nonbasis0)==1
        nonbasis0 = setdiff(nonbasis0, enteringvar);
       nonbasisu = [nonbasisu enteringvar];
    else
      nonbasisu = setdiff(nonbasisu, enteringvar);
      nonbasis0 = [nonbasis0 enteringvar];
    end
    end
        
else
basis = setdiff(basis, leavingvar);
basis = [basis enteringvar];

 if leavingbound == 0
    if ismember(enteringvar, nonbasis0)==1
    nonbasis0 = setdiff(nonbasis0, enteringvar);
    else 
    nonbasisu = setdiff(nonbasisu, enteringvar);
    end
    nonbasis0 = [nonbasis0 leavingvar];
 else
    if ismember(enteringvar, nonbasis0)==1
    nonbasis0 = setdiff(nonbasis0, enteringvar);
    else 
    nonbasisu = setdiff(nonbasisu, enteringvar);
    end
     nonbasisu = [nonbasisu leavingvar];
 end
 
end
basis = sort(basis);
 nonbasis0 = sort(nonbasis0);
 nonbasisu = sort(nonbasisu);
 B = A(:,basis);
AN0 = A(:,nonbasis0);
ANU = A(:,nonbasisu);
cb = c(basis);
cn0 = c(nonbasis0);
cnu = c(nonbasisu);
ub = u(basis);
[L,U,p] = lu(B,'vector');
y = solveyBc(L,U,p,cb);
  xu = u(nonbasisu);
if isempty(ANU) == 1
    xb = solveBxb(L,U,p,b);
else
    
    nb = b - ANU*xu;
    xb = solveBxb(L,U,p,nb);
end
rn0 = cn0 - y*AN0; 
rnu = cnu - y*ANU;
 k = k+1;
 mr0 = min(rn0);
 mru = max(rnu);
 x = zeros(n,1);
x(basis)= xb;
 xu = u(nonbasisu);
 x(nonbasisu)= xu;
z = c*x+const;
if mod(k,printlevel) == 0
   level=sprintf('Iteration =%d, z= %d',...
       k, z);
   disp(level);
end
end
end

fprintf('Phase II finished in %d iterations\n',k);
if tlimit == +Inf;
    x = zeros(n,1);
    x(basis)=xb;
   xu = u(nonbasisu);
   x(nonbasisu)= xu;
    z = -Inf;
    optbasis = basis;
    optobj = c*x;
    totaliters = k+j;
     fprintf('Detected unboundedness at phaseII iteration = %d\n',k);
else
x = zeros(n,1);
x(basis)= xb;
xu = u(nonbasisu);
x(nonbasisu)= xu;
z = c*x+const;
optobj = z;
optbasis = basis;
  totaliters = k+j;
  ray = zeros(n,1);
  test = norm(A*x-b);%check for small Ax-b
    fprintf('z=%d,Optimal solution for phaseII\n',z);
end
end

end