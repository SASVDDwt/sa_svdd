function [ P,R,f_score,temp1,TP,v2 ] = AP_SVDD_predict( label,data,apid,c_model )

%% --------------初始化---------------------------
e=0;
y=label;
x=data;
[N,C] = size(data);

FP=0;              %负类被预测成正类
TP=0;              %正类被预测成正类
TN=0;             %负类被预测成负类
FN=0;             %正类被预测成负类
temp1=[];
temp2=[];
v2=[];
%% --------------预测---------------------------
for l=1:N
    
    k=0;   %记录一个正类样本被预测成负类的次数
    m=0;  %记录一个负类样本被预测成负类的次数
    
    for i=1:apid
        
        [plabel,accuracy,decision_values]=svmpredict(y(l),x(l,:),c_model(i),'-q 0');
        
                %样本为正类
                if(y(l) == 1)
                            if(plabel~=1)
                                k=k+1;
                                
                            else%样本预测正确：即正类预测为正类
                                v2=[v2;decision_values];
                                TP=TP+1;
                                %fprintf('model:%g\n',i);
                                break;
                            end;
                 %样本为负类           
                else
                            if(plabel == -1)
                                m=m+1;
                                
                            else%样本预测错误，即负类预测为正类
                                temp2=[temp2;x(l,:)];
                                v2=[v2;decision_values];
                                FP=FP+1;
                                break;
                            end;
                end;
    end;
   
    if(k~=0)
        if(k==apid)
            FN=FN+1;
            temp1=[temp1;x(l,:)];
            v2=[v2;decision_values];
        else continue;
        end;
    else
        if(m==apid)
            TN=TN+1;
            v2=[v2;decision_values];
        else continue;
        end;
    end;
end;



%% --------------输出---------------------------


FP=FP/N;
TP=TP/N;
FN=FN/N;
TN=TN/N;

P=TP/(TP+FP);
R=TP/(TP+FN);
f_score = 2*TP/(2*TP+FP+FN);

fprintf('f-score:%g\n',f_score);
fprintf('精度:%g\n', P);
fprintf('召回率:%g\n',R);
end

