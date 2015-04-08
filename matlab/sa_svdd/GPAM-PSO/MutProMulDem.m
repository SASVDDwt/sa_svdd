function [MutNum_num,MutDe_num] = MutProMulDem(x,SelMuRat,N,D,NumOfMut,DemOfMut,interNum)
% 选择要变异的粒子及维度
%选取最聚集的维度，在选取最聚集的粒子
%%构造粒子聚集程度的矩阵D*N,clustMatrix
% clustMatrix(i,j)第i维第一个粒子在聚集范围内有多少个粒子，即多少个粒子在一堆 
% x = [1 1.9 1.4 1.8 1.3 2;3 3.1 3.1 3.8,3.7 4;5 5.5 5.5 5.9,5.8 6]';
% SelMuRat = 0.3;
% N = 6;
% D = 3;
% NumOfMut = 2;
% DemOfMut = 2;
X = x;

    % ---------------------每一维的搜索范围-------------------------------
                            DemRange = zeros(1,D);
                            clustRange = zeros(1,D);
                            for i = 1:D
                                DemRange(i) = max(X(:,i))-min(X(:,i));
                                clustRange(i) = SelMuRat*DemRange(i);
                            end 
%                             DemRange
                            clustMatrix = zeros(D,N);
                            ParMark = zeros(D,N);
                            A = cell(D,N);    % clustMatrix 每个元素对应的粒子标号
                            for dem = 1:D
                                %将每一维所有粒子按距离升序排列
                                ProDem = X(:,dem);
                                [ProDem,Index] = sort(ProDem);
                                % 重排后的粒子标号矩阵
                                ParMark(dem,:) =  Index';
                                
                                ProDemMar = Index';
                                for num = 1:N
                                    brr = find(ProDem<=ProDem(num)+clustRange(dem)) ;                                 
                                    clustMatrix(dem,num) =  length(brr)-num+1;
                                    A{dem,num} = ProDemMar(num:length(brr));                                    
                                end
                                
                            end
%                             if interNum>15 && interNum<500
%                                save(['mutaResult/AckleyFunction分布矩阵Aresult/第' num2str(interNum) '代分布矩阵Aresult.mat'],'clustMatrix'); 
%                             end
%                             if interNum>2000
%                                 disp(['第' num2str(interNum) '代聚集矩阵'])
%                                disp( clustMatrix) 
%                                disp('==================')
%                             end
%                        
                            B = zeros(1,D);     % 每一维最大堆聚集了多少个粒子
                            Bmark = cell(1,D); % 最大堆对应第几个粒子
                            for dem = 1:D
                                B(dem) = max(clustMatrix(dem,:));
                                Bmark{dem} = find(clustMatrix(dem,:) == max(clustMatrix(dem,:)));
                            end
%                              if interNum>50 && interNum<60
%                                save(['mutaResult/Penalized1分布矩阵Aresult/第' num2str(interNum) '代分布矩阵Aresult.mat'],'B'); 
%                             end
%                             B
                            %% 赌轮选择需要变异的维度
                             BMar = 1:D;
                            
                            for Mutdem = 1:NumOfMut
                                 B = B/sum(B);
                                 SelRat = cumsum(B);                         
  %  ----------------------赌轮选择得到变异的维度 DemOfMut----------------
                                SelRat = [0 SelRat];
                                a = rand;
                                for i = 1:(length(SelRat)-1)
                                     if a >SelRat(i) && a <=SelRat(i+1)
                                         Demm = i;
                                        break;
                                     end
                                end
                                MutDem_num(Mutdem) = BMar(Demm) ;
%                                       SelDem(j) =  ;
                                B(i) = 0;
                                BMar(i) = 0;
                                B = nonzeros(B)';
                                BMar = nonzeros(BMar)';
                            end
  
                            % 选择最密集堆出现的维度 即B中最大的值
%                             DemOfMut = find (B == max(B));
                            %选择对应的粒子
                            for dem = 1:NumOfMut
                            MutNum_num(dem) = A{MutDem_num(dem),Bmark{MutDem_num(dem)}(1)}(1);  
                            end
                           
                           %根据粒子找他最大的维度 
                           for Munumm = 1:NumOfMut
                               %对于每个变异粒子找维度
                               %找该粒子在ParMark中位置D个 %再找该粒子该维度上的聚集数目
                               for de = 1:D
                                   MuParMarPosit(de) = find(ParMark(de,:)==MutNum_num(Munumm)) ;
                                    MutParClustNum(de) = clustMatrix(de,MuParMarPosit(de));
                               end
                               
                               %选聚集数目最大的前DemOfMut个
                               %存放在MutDe_num
                             
                              [BBBB,index] = sort(MutParClustNum);
                              index = index';
                              MutDe_num(Munumm,:) = index(D-DemOfMut+1:D);
                              MutDe_num(Munumm,:) = fliplr(MutDe_num(Munumm,:));
                              %% 赌轮选择维度
%                             %% 维度选择： 赌轮选择
% 
%                              DNum = 1:D;
%                             
%                             for MutDem = 1:DemOfMut
%                                  MutParClustNum = MutParClustNum/sum(MutParClustNum);
%                                  SeleRat = cumsum(MutParClustNum);                         
%   %  ----------------------赌轮选择得到变异的维度 DemOfMut----------------
%                                 SeleRat = [0 SeleRat];
%                                 a = rand;
%                                 for i = 1:(length(SeleRat)-1)
%                                      if a >SeleRat(i) && a <=SeleRat(i+1)
%                                          Demm = i;
%                                         break;
%                                      end
%                                 end
%                                 MutDe_num(Munumm,MutDem) = DNum(Demm) ;
% %                                       SelDem(j) =  ;
%                                 MutParClustNum(i) = 0;
%                                 DNum(i) = 0;
%                                 MutParClustNum = nonzeros(MutParClustNum)';
%                                 DNum = nonzeros(DNum)';
%                             end
                           end
                           
end