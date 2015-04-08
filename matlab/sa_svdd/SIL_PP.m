function [pout,Sils,p2] = SIL_PP(train,s)
%{
首先对数据集进行AP聚类中P值的优化，利用SIL参数

%}


p1=max(s(:,3));
p2=median(s(:,3));
p3=min(s(:,3));

gap=(0.6*p2-8*p3)/40;
%gap=(p1-8*p3)/40;
pout=zeros(40,1);
%pout(1)=p1;
%{
for i=1:20
    pout(i)=0.6*i*p2;
end
for i=21:40
    pout(i)=0.2*i*p3;
end
%}
for i=2:40
    pout(i)=pout(i-1)-gap;
end


Sils=zeros(40,2);

for j=1:40
    [Sil]=SIL(train,s,pout(j));
    Sils(j,1)=Sil;
    
    [idx,netsim,dpsim,expref]=apcluster(s,pout(j));
    cluster=size(unique(idx),1);
    Sils(j,2) = cluster;
end

end

