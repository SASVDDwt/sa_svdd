function [ result ] = Fitness33( x,D,roted,MM,train_scale,train_label,test_scale,test_label )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 D=2;
    cmd = ['-s 2 ',' -g ',num2str( x(2) ),' -n ',num2str( x(1) )];
    model = svmtrain(train_label, train_scale, cmd);
    [ P,R,f_score,TP,TN ] = SVDD_predict( test_label,test_scale,model);
    result = f_score;
  
end