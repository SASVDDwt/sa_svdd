function [pout,p_value,Sils,cluster_num,p_percent]=plot_result( data,s )
%UNTITLED 此处显示有关此函数的摘要
%   

[pout,Sils,p2] = SIL_PP(data,s);

    p_num=find(Sils(:,1) == max(Sils(:,1)));%第几个P值
    
    cluster_num=Sils(p_num(1),2);%最佳聚类数目
       while cluster_num>20
             Sils(p_num,1) =0;
             p_num=find(Sils(:,1) == max(Sils(:,1)));
             cluster_num=Sils(p_num(1),2);%最佳聚类数目
        end
    p_value=pout(p_num(1));%对应P值
    p_percent=p_value/p2;%P值与中值关系




end

