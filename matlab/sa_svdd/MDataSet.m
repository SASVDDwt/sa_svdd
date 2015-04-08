function [ train_label,train_scale,test_label,test_scale ] = MDataSet( data,label,target )
%%
%数据处理
%输入：样本，样本标签，作为正类的标签
%随机选取正类样本的50%作为训练样本
%其余30%正类样本 + 随机负类样本的30%作为 测试样本
%输出 训练样本与测试样本


%%
%构造训练样本
ii=find(label == target);
train = data(ii,:);
[m,n] = size(train);  % m行n列矩阵（m个数据正类样本）

A = randperm(m); %产生 1~m 的一个随机排列A
M=round(m*0.7);
B = A(1:M);  %取A的前M个元素组成B，B即为M个m以内不同的随机整数；即取80%个正类样本
C = A(M+1:M+m*0.3); %取余下30%个正类样本

train_label = ones(M,1);%训练样本标签
train_scale = train(B,:); %训练样本



%%
%构造测试样本

jj = find(label ~= target);
test = data(jj,:);
[m1,n1] = size(test); % m1个负类样本

A1 = randperm(m1); %产生 1~m1 的一个随机排列A1
M1 = round(m1*0.3); 
B1 = A1(1:M1);  %取50%负类样本

test_label = [ones(round(m*0.3),1);ones(M1,1)*-1];
test_scale = [train(C,:);test(B1,:)];
end

