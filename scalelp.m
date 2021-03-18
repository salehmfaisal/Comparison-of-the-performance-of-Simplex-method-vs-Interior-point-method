%
%  [As,bs,cs,us,r,s]=scalelp(A,b,c,u)
%
% Scales a linear programming problem.
%
function [As,bs,cs,us,r,s]=scalelp(A,b,c,u)
%
% Copy the inputs to the outputs.
%
As=A;
bs=b;
cs=c;
us=u;
%
% Get the problem size.
%
[m,n]=size(As);
%
% Allocate space for the scaling vectors.
%
s=zeros(n,1);
r=zeros(m,1);
%
% Do the row scaling first.
%
for i=1:m
  r(i)=max(max(abs(As(i,:))));
  if (r(i) > 0)
    As(i,:)=As(i,:)/r(i);
    bs(i)=bs(i)/r(i);
  else
    r(i)=1.0;
  end
end
%
% Now, the column scaling.
%
for j=1:n
  s(j)=max(max(abs(As(:,j))));
  if (s(j) > 0)
    As(:,j)=As(:,j)/s(j);
    cs(j)=cs(j)/s(j);
  else
    s(j)=1;
  end
  if (u(j) < 1.0e20)
    if (s(j) > 0)
      us(j)=u(j)*s(j); 
    else
      us(j)=u(j);
    end
  else
    us(j)=u(j);
  end
end
