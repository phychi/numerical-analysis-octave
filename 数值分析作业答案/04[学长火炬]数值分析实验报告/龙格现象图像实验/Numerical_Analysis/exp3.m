clear
fa=@(x,p,q) (x.^p)*(x.^q)';
fb=@(x,f,t) f*(x.^t)';
x=[3:9];
f=[2.01 2.98 3.50 5.02 5.47 6.02 7.05];
% x=input('Please input x:');
N=length(x);
xx=x(1):0.1:x(N);
% f=input('Please input f:');
% n=input('Please input n:');
n=[1 2 4 5 6];
% a=zeros(m+1);
line={'o';'x';'+';'*';'-';':';'-.';'--';'.'};
str=cell(1,length(n));
for i=1:length(n)
    m=n(i);
    A=zeros(m+1,m+1);
    b=zeros(m+1,1);
    for j=1:m+1
        b(j)=fb(x,f,j-1);
        for k=1:m+1
            A(j,k)=fa(x,j-1,k-1);
        end
    end
    a=inv(A)*b;
    y=zeros(1,length(xx));
    for j=1:m+1
        y=y+a(j)*(xx.^(j-1));
    end
    ee=y([1:10:61])-f;
    e=ee*ee'
    plot(xx,y,line{randi(length(line))},'Color',rand(1,3));
    str{i}=[num2str(m) ' polynomial fitting curve'];
    hold on;
end
str{length(n)+1}='raw';
scatter(x,f,'r','filled');
legend(str,'FontSize',16,'Location','NorthWest');
ylim([2 7.5]);
grid on;