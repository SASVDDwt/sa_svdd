function  [f_score1,f_score2,f_score3,f_score4]=mainSelfAdaptiveAPSVDD(filename,train_label,train_scale,test_label,test_scale,s,ppp)
demantion = 2;
accuracy = 10^(-5);
c3 = 0.4;
w = 0.4;
SelMuRat = 0.1;
roted = 0;


fid=fopen(filename,'a');
fprintf(fid,'------------------------------\n');

%%  svdd+ap
pp=zeros(10,1);
rr=zeros(10,1);
ff=zeros(10,1); 
aa=zeros(10,1);
fprintf(fid,'SVDD+AP\n');

%[s,p] = SIL_P(train_scale);
for i=1:10
 [apid,c_model] = AP_SVDD_train2( train_label,train_scale,test_label,test_scale,s,ppp);
 [ P,R,f_score,temp1,TP] = AP_SVDD_predict( test_label,test_scale,apid,c_model );
 pp(i)=P;
 rr(i)=R;
 ff(i)=f_score;
 aa(i)=apid;
 
 fprintf(fid,'-%d',i);
 fprintf(fid,'-\n');
 fprintf(fid,'P:%f\n', P);
 fprintf(fid,'R:%f\n', R);
 fprintf(fid,'f_score:%f\n\n', f_score);
 
end
 fprintf(fid,'-平均值-\n');
 fprintf(fid,'mean-P:%f\n', mean(pp));
 fprintf(fid,'mean-R:%f\n', mean(rr));
 fprintf(fid,'mean-F:%f\n\n', mean(ff));
f_score1=mean(ff);

%% svdd

pp=zeros(10,1);
rr=zeros(10,1);
ff=zeros(10,1); 
fprintf(fid,'~~~~~~~~~~~~~~~~~SVDD\n');
for i=1:10
resutl = PSOLHSimproveZhongMuMinQHMinDmF(@Fitness11,@regress11,@LHspso,10,0.9,0.9,c3,w,80,2,'Penalized1',50,accuracy,SelMuRat,roted,train_scale,train_label,train_scale,train_label);
option= ['-s 5',' -g ',num2str(resutl.xm(2)),' -c ',num2str(resutl.xm(1))];
%[bestCVaccuarcy,bestc,bestg,pso_option] = psoSVMcgForClassCV(train_label,train_scale);
%option= ['-s 5',' -g ',num2str(bestg),' -c ',num2str(bestc)];

model=svmtrain(train_label,train_scale,option);
 [ P,R,f_score,TP,TN ] = SVDD_predict( test_label,test_scale,model);
 pp(i)=P;
 rr(i)=R;
 ff(i)=f_score;

 fprintf(fid,'-%d',i);
 fprintf(fid,'-\n');
 fprintf(fid,'P:%f\n', P);
 fprintf(fid,'R:%f\n', R);
 fprintf(fid,'f_score:%f\n\n', f_score);
 
end
 fprintf(fid,'-平均值-\n');
 fprintf(fid,'mean-P:%f\n', mean(pp));
 fprintf(fid,'mean-R:%f\n', mean(rr));
 fprintf(fid,'mean-F:%f\n\n', mean(ff));
 f_score2=mean(ff);

%%  svdd+kmeans
pp=zeros(10,1);
rr=zeros(10,1);
ff=zeros(10,1); 
fprintf(fid,'~~~~~~~~~~~~~SVDD+KMEANS\n');
for i=1:10
 [ apid,c_model ] = Kmeans_SVDD( train_label,train_scale,mean(aa),test_label,test_scale)
[ P,R,f_score,temp1,TP]  = AP_SVDD_predict( test_label,test_scale,apid,c_model );

pp(i)=P;
 rr(i)=R;
 ff(i)=f_score;

 fprintf(fid,'-%d',i);
 fprintf(fid,'-\n');
 fprintf(fid,'P:%f\n', P);
 fprintf(fid,'R:%f\n', R);
 fprintf(fid,'f_score:%f\n\n', f_score);
end
 fprintf(fid,'-平均值-\n');
 fprintf(fid,'mean-P:%f\n', mean(pp));
 fprintf(fid,'mean-R:%f\n', mean(rr));
 fprintf(fid,'mean-F:%f\n\n', mean(ff));
 f_score3=mean(ff);
 
%% one-class

pp=zeros(10,1);
rr=zeros(10,1);
ff=zeros(10,1); 
fprintf(fid,'~~~~~~~~~~~~~~~ONE-CLASS\n');
for i=1:10
%需修改psoSVMcgForClass中参数
resutl = PSOLHSimproveZhongMuMinQHMinDmF(@Fitness22,@regress11,@LHspso,10,0.9,0.9,c3,w,80,2,'Penalized1',50,accuracy,SelMuRat,roted,train_scale,train_label,train_scale,train_label);
option= ['-s 2',' -g ',num2str(resutl.xm(2)),' -n ',num2str(resutl.xm(1))];
%[bestCVaccuarcy,bestc,bestg,pso_option] = psoSVMgnForClassCV(train_label,train_scale);
%option= ['-s 2',' -g ',num2str(bestg),' -n ',num2str(bestc)];

model=svmtrain(train_label,train_scale,option);
 [ P,R,f_score,TP,TN ] = SVDD_predict( test_label,test_scale,model);
pp(i)=P;
 rr(i)=R;
 ff(i)=f_score;

 fprintf(fid,'-%d',i);
 fprintf(fid,'-\n');
 fprintf(fid,'P:%f\n', P);
 fprintf(fid,'R:%f\n', R);
 fprintf(fid,'f_score:%f\n\n', f_score);
end
 fprintf(fid,'-平均值-\n');
 fprintf(fid,'mean-P:%f\n', mean(pp));
 fprintf(fid,'mean-R:%f\n', mean(rr));
 fprintf(fid,'mean-F:%f\n\n', mean(ff));
 f_score4=mean(ff);
  fid=fclose(fid);
end

