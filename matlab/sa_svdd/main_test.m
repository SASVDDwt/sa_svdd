function [ f_score ] = main_test( train_label,train_scale,test_label,test_scale )


s=s_value(train_scale);
[pout,p_value,Sils,cluster_num,p_percent]=plot_result( train_scale,s );


ff=zeros(10,1); 
for i=1:10
 [apid,c_model] = AP_SVDD_train2( train_label,train_scale,test_label,test_scale,s,p_value);
 [ P,R,f_score,temp1,TP] = AP_SVDD_predict( test_label,test_scale,apid,c_model );

 ff(i)=f_score;
end


f_score = mean(ff);
end

