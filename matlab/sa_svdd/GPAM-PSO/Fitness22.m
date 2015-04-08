function [ result ] = Fitness22( x,D,roted,MM,train_scale,train_label,test_scale,test_label )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 D=2;
 
    cmd = ['-s 2 ',' -v 3',' -g ',num2str( x(2) ),' -n ',num2str( x(1) )];
   
    result = svmtrain(train_label, train_scale, cmd);
  
end