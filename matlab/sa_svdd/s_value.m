function  s=s_value( data )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% ---------------计算数据集的相似度值----------------------

x=data;

[N,C]=size(x);
%[N1,C1]=size(test_scale);
% ---构造相似度矩阵---  

M=N*N-N;
s=zeros(M,3);
j=1;
for i=1:N
  for k=[1:i-1,i+1:N] %k 的取值是1~N里除了 i的所有整数
    s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);
    j=j+1;
  end;
end; 


% ---分别求S的最大值、最小值、中位数、四分位数、第三四分位数----

end
