function [ apid,c_model ] = Kmeans_SVDD( label,data,k,test_label,test_scale)
%--------------初始化------------------  
demantion = 2;
accuracy = 10^(-5);
c3 = 0.4;
w = 0.4;
SelMuRat = 0.1;
roted = 0;

t1=clock;
y=label;
x=data;
[N1,M1]=size(test_scale);
%--------------kmeans 聚类------------------  
[idx]=kmeans(data,k);


%--------------初始化聚类个数，支持向量数------------------  

apid=0;
SVS=0;
%sv=[];

%--------------对每个聚类结果进行SVDD训练-----------------

for i=unique(idx)'
  ii=find(idx==i);  
  apid=apid+1;
  
    %--------------PSO为每个子类进行参数寻优-----------------  
%[bestCVaccuarcy,bestc,bestg,pso_option] = psoAPSVMcg(y(ii,1),x(ii,:),test_label,test_scale,test_label1,test_scale1);
% [bestCVaccuarcy,bestc,bestg,pso_option] = psoSVMcgForClass(y(ii,1),x(ii,:),test_label,test_scale);
%[bestacc,bestc,bestg] = SVMcgForClass(y(ii,1),x(ii,:));
  %--------------END PSO-----------------  
   resutl = PSOLHSimproveZhongMuMinQHMinDmF(@Fitness11,@regress11,@LHspso,10,0.9,0.9,c3,w,80,2,'Penalized1',50,accuracy,SelMuRat,roted,x(ii,:),y(ii,1),x(ii,:),y(ii,1));
  option= ['-s 5',' -g ',num2str(resutl.xm(2)),' -c ',num2str(resutl.xm(1))];
  %option= ['-s 5',' -g ',num2str(bestg),' -c ',num2str(bestc)];
  %SVDD
   %option=['-s 5 -g 0.7755   -c 0.1818'];
  model=svmtrain(y(ii,1),x(ii,:),option);
  c_model(apid)=model;
  SVS=SVS+model.totalSV;
  % sv=[sv;model.SVs]; 
  
end;


%--------------输出----------------- 


fprintf('kmeans result:%d\n',apid);
fprintf('total SVs:%d\n',SVS);

t2=clock;
fprintf('run time:%g\n',etime(t2,t1));


 %--------------PSO为每个子类进行参数寻优-----------------  
 %[bestCVaccuarcy,bestc,bestg,pso_option] = psoAPSVMcg(ones(SVS,1),sv,test_label,test_scale,test_label,test_scale);
  %--------------END PSO-----------------  
  
 % option= ['-s 5',' -c ',num2str(bestc ),' -g ',num2str(bestg )];
 % model=svmtrain(ones(SVS,1),sv,option);
%model=svmtrain(ones(SVS,1),sv,'-s 5');
%[plabel,accuracy]=svmpredict(label,data,model);
end
