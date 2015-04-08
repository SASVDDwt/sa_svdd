function [ Sil ] = SIL( data,s,p )
%{
    1.ap聚类
    2.构造silhouette中clust参数，为每个数据所属类别的矩阵
                  根据ap聚类得到的idx矩阵，为所有聚类中心重新赋值（赋值为1到所有聚类个数）
    3.计算silhouette值： R = silhouette(data,classlabel,'Euclidean'); 
                  data:数据集
                  classlabel:每个样本点所属类别的矩阵
                  'Euclidean':欧式距离
                   R:返回值
    4.计算所有样本点SIL值的均值
%}
[idx,netsim,dpsim,expref]=apcluster(s,p);
classlabel=idx;
cluster=1;
for i=unique(idx)'
    cluster_label= i == idx;
    classlabel(cluster_label) = cluster;
    cluster=cluster+1;
end
 R = silhouette(data,classlabel,'Euclidean');
 Sil=mean(R);
end

