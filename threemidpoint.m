function [inthree,iterthree] =threemidpoint(a,b,m,fun,tol,alpha,int)
inthree=0;
iterthree=0;
while abs(inthree-exp(3)+1)>=tol
h=(b-a)/m;
x=a+h/2:h:b;
dim=length(x);
y=feval(fun,x);
if size(y)==1
    y=diag(ones(dim))*y;
end
inthree=1/3*(int+h*sum(feval(fun,alpha+h+h/6))+feval(fun,alpha+h+5*h/6));
m=m*2;
iterthree=iterthree+1;
if abs(inthree-exp(3)+1)<=tol
    break;
end
end
return 