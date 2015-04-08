function [f_score1,f_score2,f_score3,f_score4]=main(input_filename,target,output_filename)
%MAIN Summary of this function goes here
[label,data]=libsvmread(input_filename);
[ train_label,train_scale,test_label,test_scale ] = MDataSet( data,label,target );
s=s_value(train_scale);
[pout,p_value,Sils,cluster_num,p_percent]=plot_result( train_scale,s );

  [f_score1,f_score2,f_score3,f_score4]=mainSelfAdaptiveAPSVDD(output_filename,train_label,train_scale,test_label,test_scale,s,p_value);
end

