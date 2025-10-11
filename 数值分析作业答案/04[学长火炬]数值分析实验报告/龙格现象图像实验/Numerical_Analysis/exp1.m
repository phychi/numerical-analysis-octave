clear all
%% 二分法
f=@(x) sin(x)-x^2/2;
r=[1 2];
root1=0;
k=0;
while r(2)-r(1)>=(0.5e-5)
    rm=mean(r);
    k=k+1;
    if f(r(1))*f(rm)<0
        r(2)=rm;
    end
    if ~ (f(r(1))*f(rm))
            root=rm;
            break;
    end
    if f(r(1))*f(rm)>0
        r(1)=rm;
    end
    root1=mean(r);
end
root1
k
f(root1)
%% 牛顿法1
fn=@(x) x-(x*exp(x)-1)/(x*exp(x)+exp(x));
x0=0.5;
x1=fn(x0);
k=1;
while abs(x1-x0)>=(0.5e-5)
    k=k+1;
    x0=x1;
    x1=fn(x0);
end
x1
k
x1*exp(x1)-1
%% 牛顿法2
f=@(x) x-(x^3-x-1)/(3*x^2-1)
x0=1;
x1=f(x0);
k=1;
while abs(x1-x0)>=(0.5e-5)
    x0=x1;
    x1=f(x0);
    k=k+1;
end
x1
k
x1^3-x1-1
%% 牛顿法3_1
f=@(x) x-((x-1).^2.*(2*x-1))./(2*(x-1).*(2*x-1)+(x-1).^2*2);
x0=0.45;
x1=f(x0);
k=1;
while abs(x1-x0)>=(0.5e-5)
    x0=x1;
    x1=f(x0);
    k=k+1;
end
x1
k
(x1-1)^2*(2*x1-1)
%% 牛顿法3_2
f=@(x) x-((x-1).^2.*(2*x-1))./(2*(x-1).*(2*x-1)+(x-1).^2*2);
x0=0.65;
x1=f(x0);
k=1;
while abs(x1-x0)>=(0.5e-5)
    x0=x1;
    x1=f(x0);
    k=k+1;
end
x1
k
(x1-1)^2*(2*x1-1)
%% 割线法
g=@(x) x*exp(x)-1;
f=@(x0,x1) x1-g(x1)/(g(x1)-g(x0))*(x1-x0)
r=[0.4 0.6];
k=0;
while abs(r(1)-r(2))>=(0.5e-5)
    x=f(r(1),r(2));
    r(1)=r(2);
    r(2)=x;
    k=k+1;
end
r(2)
k
x*exp(x)-1
%% 改进牛顿法
fn=@(x) x-2*((x-1).^2.*(2*x-1))./(2*(x-1).*(2*x-1)+(x-1).^2*2);
x0=0.55;
x1=fn(x0);
k=1;
while abs((x1-1)^2*(2*x1-1))>=(0.5e-5)  %若取abs(x1-x0)>=0.5e-5，无法满足
    x0=x1;
    x1=fn(x0);
    k=k+1;
end
x1
k
(x1-1)^2*(2*x1-1)
%% 拟牛顿法
A=@(x,y,z) inv([y x -2*z
    y*z-2*x x*z+2*y x*y 
    exp(x) -exp(y) 1]);
F=@(x,y,z) [x*y-z^2-1
    x*y*z+y^2-x^2-2
    exp(x)+z-exp(y)-3];
f=@(x,y,z) [x;y;z]-A(x,y,z)*F(x,y,z);
x=1;
y=1;
z=1;
a=f(x,y,z);
k=1;
while norm(F(a(1),a(2),a(3)))>=0.5e-5
    a=f(x,y,z);
    x=a(1);
    y=a(2);
    z=a(3);
    k=k+1;
end
[x y z]
k
norm(F(a(1),a(2),a(3)))