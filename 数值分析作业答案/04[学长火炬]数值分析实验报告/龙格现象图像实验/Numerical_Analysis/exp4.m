clear
f=@(x) sin(x.^2);
% mm=input('请输入积分上下限：');
mm=[0,1];
e=input('请输入精度：');
H=mm(2)-mm(1);
T(1,1)=H*(f(mm(1))+f(mm(2)))/2;
x=mean(mm);
T(2,1)=T(1,1)+H*f(x)/2;
T(2,2)=(4*(T(2,1))-T(1,1))/(4-1);
i=2;
while abs(T(size(T,1),size(T,2))-T(size(T,1),size(T,2)-1)) >= e
    i=i+1;
    T(i,1)=T(i-1,1)/2+H/(2^(i-1))*sum(f(mm(1)+H/(2^(i-1)):H/(2^(i-2)):mm(2)));
    for j=2:i
        T(i,j)=((4^(i-1))*T(i,j-1)-T(i-1,j-1))/(4^(i-1)-1);
    end
end
T
abs(T(size(T,1),size(T,2))-T(size(T,1),size(T,2)-1))
T(size(T,1),size(T,2))