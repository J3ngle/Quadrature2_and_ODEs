%% FCM Project 4
%Jonathan Engle
%Due March 12
%% Define Functions and parameters
clear;clc;clear all 
tic
format short
%lower bound 
global a
a=0;
%upper bound
global b
b=3;
% Global alpha
global alpha
alpha=1/3;
%Global tol
global tol 
%tol=0.01; %0.0001
tol=.0001;

%function
global fun 
fun=@(t)(exp(t));
%% Composite trapezpodal rule
[inttrap,mtrap] =trap(a,b,fun,tol);
error=abs(inttrap-exp(b)+1);
mtrap;
%% Global refinement CTR
[twoint,prev_integral,n,eval] = ctr(fun, a, b, tol);
twoint;
twoint-exp(3)+1;
n;
eval;
r=log(((twoint-prev_integral)/(exp(3)-1)-twoint)+1)/log(2);
%% Composite Gauss legendre method
% x=@(z)exp((z*(b-a)+(b+a)/2));
% dx=(b-a)/2;
z0=-1/sqrt(3);
z1=1/sqrt(3);

m=1;
hm=(b-a)/m;
ai=linspace(a,b-hm,m);
% x=@(z)((b-a)/2*z+(b+a)/2);
% F=@(t)fun(x(t))*(b-a)/2;
% Gauss=F(z0)+F(z1);
Gauss=0;
while abs(Gauss-exp(3)+1)>=tol
    Gauss=(hm)/2*(exp((hm)/2*(z0)+hm/2+ai(1))+exp((hm)/2*(z1)+(hm)/2+ai(1)));
    for i=2:m
        Gauss = Gauss+(hm)/2*(exp((hm)/2*z0+hm/2+ai(i))+exp((hm)/2*(z1)+(hm)/2+ai(i)));
    end
    m=m+1;
    hm=(b-a)/m;
    ai=linspace(a,b-hm,m);
end
m;
%% Forward and Backward Euler
%Set up
lambda=1; % 1,-1,-0.01
omega=0.0001; %10,0.01
eulerfun = @(t, y)lambda*y;% lambda*(y-sin(omega*t))+omega*cos(omega*t);%lambda*t; % %Manually enter functions

%IC's that are indexed +1
t(1) = 0;
y(1) = 1;
%Final time
T = 10; 
%Step sze
h = .0001;      
%itterations
N = ceil((T - t(1)) / h);
t = zeros(N+1, 1);
y = zeros(N+1, 1);
t(1) = 0;
y(1) = 1;
true_eulerfun=@(t,y) exp(t*lambda);%@(t,y) @(t,y) (y(1))*(sin(0))*exp(lambda*t) +sin(omega*t);%exp(t*lambda);%

% Forward Euler method
for i = 1:N
    t(i+1) = t(i) +h;     
    y(i+1) = y(i) + h * eulerfun(t(i), y(i));
end
%Change t and y for backward
bt = zeros(N+1, 1);
by = zeros(N+1, 1);
bt(1) = 0;
by(1) = 1;
% Backward Euler method
for i = 1:N
    bt(i+1) = bt(i) + h;     
    by(i+1) = by(i) + h * eulerfun(t(i+1), y(i+1));
end

g=linspace(t(1),10);
et=exp(g);
norm(et-y,inf);
Error_Forward=norm(feval(true_eulerfun,t,y)-y,inf)
Error_Backward=norm(feval(true_eulerfun,bt,by)-by,inf)
hold on
plot(t,y,'b->')
plot(bt,by,'m-<')
%ode45(true_eulerfun,0:10,0)
%plot(exp(-t))
legend({'Forward Euler','Backward Euler'},'Location','northeastoutside')
xlabel('time')
ylabel('y-axis')
hold off
toc

