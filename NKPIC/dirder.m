function z = dirder(x,w,f,f0)
% Finite difference directional derivative
% Approximate f'(x) w
% 
% C. T. Kelley, November 25, 1993
%
% This code comes with no guarantee or warranty of any kind.
%
% function z = dirder(x,w,f,f0)
%
% inputs:
%           x, w = point and direction
%           f = function
%           f0 = f(x), in nonlinear iterations
%                f(x) has usually been computed
%                before the call to dirder

%
% Hardwired difference increment.
epsnew=1.d-7;
%
n=length(x);
%
% scale the step
%
if norm(w) == 0
    z=zeros(n,1);
return
end
epsnew = epsnew/norm(w);
if norm(x) > 0
    epsnew=epsnew*norm(x);
end
%
% del and f1 could share the same space if storage
% is more important than clarity
%
del=x+epsnew*w;
f1=feval(f,del);
z = (f1 - f0)/epsnew;
