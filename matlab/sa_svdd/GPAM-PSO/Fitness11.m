function [ result ] = Fitness11( x,D,roted,MM,train_scale,train_label,test_scale,test_label )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
 D=2;
 
    cmd = ['-s 5 ',' -v 3',' -g ',num2str( x(2) ),' -c ',num2str( x(1) )];
   
    result = svmtrain(train_label, train_scale, cmd);
  
end

