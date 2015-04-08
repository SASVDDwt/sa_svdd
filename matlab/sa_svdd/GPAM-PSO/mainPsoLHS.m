function mainPsoLHS
% demList = [200 500 1000];
% for dem = 1:3
%     demantion = demList(dem);


%%demantion = 30;
demantion = 2;
accuracy = 10^(-5);
c3 = 0.4;
w = 0.4;
SelMuRat = 0.1;
roted = 0;
%     for SelMuRat = 0.3:0.1:0.7
%     for w = 0.2:0.2:0.8;
%     for c3 = 0:0.2:2
% -----------all the test function ---------------------------------------------------------------------------------------------------------------------------
% % %  ----------------µ•∑Â---------------------------------------------
%         PSOLHSimproveZhongMuMinQHMinDmF(@Sphere,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Sphere',100,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Schwefel,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Schwefel',10,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Quardirc,@regress11,@LHspso,30,0.9,0.9,c3,w,20000,demantion,'Quardirc',100,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Schwefel2,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Schwefel2',100,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Rosenbrock,@regress11,@LHspso,30,0.9,0.9,c3,w,2000,demantion,'Rosenbrock',2,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Step,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Step',100,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Tablet,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Tablet',100,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Ellipse,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Ellipse',100,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@DiffPower,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'DiffPower',100,accuracy,SelMuRat,roted);
% % % ------------------∂‡∑Â-----------------------------------------------
%         PSOLHSimproveZhongMuMinQHMinDmF(@SchwefelMaxFunction,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'SchwefelMaxFunction',500,accuracy,SelMuRat,roted);
% PSOLHSimproveZhongMuMinQHMinDmF(@Dminima,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Dminima',5,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Rastrigin,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Rastrigin',5,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@WeiMaxFunction,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'WeiMaxFunction',5.12,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@AckleyFunction,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'AckleyFunction',32,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Griewank,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Griewank',600,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@NonMaxFunction,@regress11,@LHspso,30,0.9,0.4,c3,w,5000,demantion,'NonMaxFunction',0.5,accuracy,SelMuRat,roted);
% %         PSOLHSimproveZhongMuMinQHMinDmF(@Salomon,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Salomon',100,accuracy,SelMuRat,roted);
% %          PSOLHSimproveZhongMuMinQHMinDmF(@SchwefelInteger,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'SchwefelInteger',pi,accuracy,SelMuRat,roted);
%%result = PSOLHSimproveZhongMuMinQHMinDmF(fitness,regress,LHS,N,c1,c2,c3,w,M,D,fileName,boundary,accuracy,SelMuRat,roted,train_scale,test_scale)
         PSOLHSimproveZhongMuMinQHMinDmF(@Fitness,@regress11,@LHspso,10,0.9,0.9,c3,w,200,2,'Penalized1',50,accuracy,SelMuRat,roted,train_scale,test_scale);
         PSOLHSimproveZhongMuMinQHMinDmF(@Fitness,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Penalized2',50,accuracy,SelMuRat,roted,train_scale,test_scale);
% % -------------------‘Î…˘---------------------------------------------------
%         PSOLHSimproveZhongMuMinQHMinDmF(@NoisySchwefel,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'NoisySchwefel',100,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@NoisyQuartic,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'NoisyQuartic',100,accuracy,SelMuRat,roted);
% -------------------mis-scaled---------------------------------------------
%         PSOLHSimproveZhongMuMinQHMinDmF(@Rastrigin10,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Rastrigin10',5,accuracy,SelMuRat,roted);
%         PSOLHSimproveZhongMuMinQHMinDmF(@Rastrigin100,@regress11,@LHspso,30,0.9,0.9,c3,w,5000,demantion,'Rastrigin100',5,accuracy,SelMuRat,roted);

%     end
end
% end