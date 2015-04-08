function [ result ] = Fitness( x,D,roted,MM,train_scale,train_label,test_scale,test_label )
   
    D=2;
    cmd = ['-s 5 ',' -g ',num2str( x(2) ),' -c ',num2str( x(1) )];
    model = svmtrain(train_label, train_scale, cmd);
    [ P,R,f_score,TP,TN ] = SVDD_predict( test_label,test_scale,model);
    result = f_score;
   

end

