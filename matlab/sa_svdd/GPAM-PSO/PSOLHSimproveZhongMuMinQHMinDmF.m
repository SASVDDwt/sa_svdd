function result = PSOLHSimproveZhongMuMinQHMinDmF(fitness,regress,LHS,N,c1,c2,c3,w,M,D,fileName,boundary,accuracy,SelMuRat,roted,train_scale,train_label,test_scale,test_label)

%  最后对于Dminima比较好的,
% vmax = 0.9;c1=c2=c3 = 0.9;w = 0.7
% 变异策略为胡老师，< 目前均值加本身变异粒子适应度
% 不好则其他变异。
% 修改变异策略，怎样选择粒子，怎样选择维度。
% 100---4000加入minwhere。变异初具规模后 隔几十代变异
% 赌轮选择要变异的某些粒子的某些维 PSOLHS2；连续失败10次 变异变异3次选最好的 前100代 minwhere，后100代 minwhere
% N:粒子数量；
% M:迭代代数；
% D:粒子维数；
% Q:运行次数；
% for
firstNumMinwhere = 50;                     % 前多少代加minwhere；
eachIntMinwhere = 100;                     % 多少代加minwhere
lastIntMinwhere = 50;                       % 最后多少代加minwhere
%---------------参数初始化-----------------
% min_dis = 0.05;
% accuMark = 0;
% for rand11 = 0.01:0.1:1


%---------------c，g参数范围-----------------
cmin=0.001;
cmax=1;
gmin=0.001;
gmax=10;
boundary=[cmin cmax;gmin gmax];
%---------------end-----------------

Q=1;
failSeqNumMax = 5;                          % 连续失败上限
failSeqNum = 0;                             % 无显著提高则累加失败次数
num = 0;
cn1 = 0.2;
cn2 = 0.2;
num1 = 0;                                   % 一共变异的次数（每次有多个粒子）
% num11 = 0;                                  % 提高了的粒子数
% num1111 = 0;
numNormRandMin_where = 4;                   %  min_where正态分布多少点
p=zeros(1,N);                               %  历史最优适应度值p(i)
y=zeros(N,D);                               %  历史最优适位置  y(i)
Pbest = zeros(1,M);                         %  某代全局最优适应度值
min_whereFitness = zeros(1,M);
vmax = 0.9;
result.Pbest = [];
result.runtime = [];
result.xm = [];
result.fv = [];
result.accept_iter = [];
result.successEvaluation = [];              %  达到精度的适应度值计算数
result.figPbest = [];                       %  Q次最好结果的Pbest 画图
result.min_whereFitness = [];
regressX = [];
regressXpca = [];
result.x = [];
result.fitness = [];
result.min_wherewhere = [];
intervalNum = 30;
% clusMatrix = zeros(D,intervalNum);
% SelMuRat = 0.3;
% suus = 0;
%--------------------
DemOfMut = 2;
NumOfMut = 2;
%--------------------
rand11 = 0.7;
rand22 = 0.2;
format long;
MM = gallery('randhess',D,'double');

%% Q次运行程序 Q次 /

disp(fileName);
result.firstNumPbest = [];                                         %  记录前n代加入minwhere后的最后pbest

for runNum = 1:Q
    pcaMatrix = [];
    MuDistance = [];
    %             eee = 0;
    wheel = 0;
    accept_iter = 0;                                                %  达到精度的代数
    sumEvaluation = 0;
    successEvaluation = 0;
    NormRnd = zeros(numNormRandMin_where,D);                        %  min_where 正态分布点坐标
    
    %%%%%%%%%%%%%%%%%%为了更新公式方便
    MutNum_num = ones(1,NumOfMut)*400;
    MutEither = 0;
    %%%%%%%%%%%%%%%
    disp([ num2str(runNum) '=====================']);
    
    %-----------random creat,based on the clock------------------
    stream = RandStream('mt19937ar', 'Seed', sum(100 * clock));
    RandStream.setGlobalStream(stream);
    tic;                                                           % tic计时
    
    %% 初始化种群的个体位置和速度
    %             v = rand(N,D) * (2*boundary) - boundary;
    v = rand(N,D);
    x = LHS(N,D);
    % x = (2*boundary)*rand(N,D)-boundary;
    
    %% 初始化 局部最优 全局最优
    for i=1:N
        p(i)=fitness(x(i,:),D,roted,MM,train_scale,train_label,test_scale,test_label);                                    % p(i)为历史最优适应度值
        y(i,:)=x(i,:);                                             % y(i)为历史最优位置
    end
    pg = x(N,:);
   disp([pg '1=====================']);
    for i=1:(N-1)
        x1=fitness(x(i,:),D,roted,MM,train_scale,train_label,test_scale,test_label);
        x2=fitness(pg,D,roted,MM,train_scale,train_label,test_scale,test_label);
        
        if fitness(x(i,:),D,roted,MM,train_scale,train_label,test_scale,test_label)>fitness(pg,D,roted,MM,train_scale,train_label,test_scale,test_label)
            pg=x(i,:);                                             % Pg 为全局最优位置
        end
    end
    disp([pg '2=====================']);
    %min_where = rand(1,D)*(2*boundary) - boundary;                % 初始化 min_where 位置
    min_where = ([(cmax-cmin) (gmax-gmin)]) - [(rand*cmax-cmin) rand*(gmax-gmin)];
    %% 按照公式迭代 M代
    for interNum=1:M
        clusMatrix = zeros(D,intervalNum);
        %                 SelDem = [];           % 变异的维度
        % =============================更新N个粒子================================= ;
        for MarkOfPar=1:N
            
            if interNum >firstNumMinwhere && interNum<M-lastIntMinwhere
                if mod(interNum,eachIntMinwhere) == 0
                    
                    v(MarkOfPar,:)=w*v(MarkOfPar,:)+c1*rand*(y(MarkOfPar,:)-x(MarkOfPar,:))...
                        +c2*rand*(pg-x(MarkOfPar,:))+c3*rand*(min_where-x(MarkOfPar,:));
                    x(MarkOfPar,:)=x(MarkOfPar,:)+v(MarkOfPar,:);
                else
                    if MutEither == 1 &&  ~all( MutNum_num-MarkOfPar)
                        %     disp(['&&&&&&&&&&&&' num2str(MarkOfPar)])
                        v(MarkOfPar,:)=w*v(MarkOfPar,:)+cn1*rand*(y(MarkOfPar,:)...
                            -x(MarkOfPar,:))+cn2*rand*(pg-x(MarkOfPar,:));
                        x(MarkOfPar,:)=x(MarkOfPar,:)+v(MarkOfPar,:);
                    else
                        v(MarkOfPar,:)=w*v(MarkOfPar,:)+c1*rand*(y(MarkOfPar,:)...
                            -x(MarkOfPar,:))+c2*rand*(pg-x(MarkOfPar,:));
                        x(MarkOfPar,:)=x(MarkOfPar,:)+v(MarkOfPar,:);
                        
                    end
                    %
                    %                             if abs(interNum-100*(ceil(interNum/100)))<5
                    %                               disp(['第' num2str(interNum) '代' '第' num2str(MarkOfPar) '个粒子'])%
                    %                               disp('速度' )
                    %                               v(MarkOfPar,:)
                    %                               disp('位置')
                    %                               x(MarkOfPar,:)%
                    %                             end
                    
                end
            elseif  interNum<=firstNumMinwhere
%size (v)
%size (x)
%size (y)
%size (pg)
%size (min_where)

                
                v(MarkOfPar,:)=w*v(MarkOfPar,:)+c1*rand*(y(MarkOfPar,:)-x(MarkOfPar,:))...
                    +c2*rand*(pg-x(MarkOfPar,:))+c3*rand*(min_where-x(MarkOfPar,:));
                x(MarkOfPar,:)=x(MarkOfPar,:)+v(MarkOfPar,:);
                %                         elseif interNum>=4000
                %                             v(MarkOfPar,:)=w*v(MarkOfPar,:)+c2*rand*(pg-x(MarkOfPar,:))+c3*rand*(min_where-x(MarkOfPar,:));
                %                             x(MarkOfPar,:)=x(MarkOfPar,:)+v(MarkOfPar,:);
            elseif interNum>=M-lastIntMinwhere
                
                if fitness(min_where,D,roted,MM,train_scale,train_label,test_scale,test_label)<fitness(pg,D,roted,MM,train_scale,train_label,test_scale,test_label)
                    %                             disp('===========')+c2*rand*(pg-x(MarkOfPar,:))
                    v(MarkOfPar,:)=w*v(MarkOfPar,:)+c1*rand*(y(MarkOfPar,:)-x(MarkOfPar,:))...
                        +c3*rand*(min_where-x(MarkOfPar,:));
                    x(MarkOfPar,:)=x(MarkOfPar,:)+v(MarkOfPar,:);
                else
                    %                                 disp('_____________')+c3*rand*(pg-x(MarkOfPar,:))
                    v(MarkOfPar,:)=w*v(MarkOfPar,:)+c1*rand*(y(MarkOfPar,:)-x(MarkOfPar,:))...
                        +c2*rand*(pg-x(MarkOfPar,:));
                    x(MarkOfPar,:)=x(MarkOfPar,:)+v(MarkOfPar,:);
                end
            end
            
            
            %% ----------------------------边界鉴定-----------------------------------
            if x(MarkOfPar,1)<cmin;
                    x(MarkOfPar,1)=cmin+rand;
                    num = num +1;
                    disp('1==========')
            end
            if x(MarkOfPar,1)>cmax;
                    x(MarkOfPar,1)=rand*cmax;
                    num = num +1;
                     disp('2==========')
            end
             if x(MarkOfPar,2)<gmin;
                    x(MarkOfPar,2)=gmin+rand;
                    num = num +1;
                     disp('3==========')
            end
            if x(MarkOfPar,2)>gmax;
                    x(MarkOfPar,2)=rand*gmax;
                    num = num +1;
                     disp('4==========')
            end
            x(MarkOfPar,2)
            x(MarkOfPar,1)
           %{
                for k=1:D
                if x(MarkOfPar,k)<-boundary;
                    x(MarkOfPar,k)=-boundary*rand;
                    %                            if interNum>50
                    %                            disp(['第几个粒子' num2str(MarkOfPar) '第几维' num2str(k)])
                    %                            disp('==========')
                    %                            disp(interNum);
                    %                            end
                    num = num +1;
                elseif x(MarkOfPar,k)>boundary
                    x(MarkOfPar,k)=rand*boundary;
                    num = num +1;
                end
            end
            %}
            %------------------------------速度鉴定----------------------------------
            for k = 1:D
                if v(MarkOfPar,k)>vmax
                    %                            num = num +1;
                    v(MarkOfPar,k)=rand*vmax;
                elseif v(MarkOfPar,k)<-vmax
                    %                             num = num +1;
                    v(MarkOfPar,k)=-vmax*rand;
                end
            end
            %                       if abs(interNum-100*(ceil(interNum/100)))<5
            %                               disp(['第' num2str(interNum) '代验证后' '第' num2str(MarkOfPar) '个粒子'])
            %                               v(MarkOfPar,:)
            %                               x(MarkOfPar,:)
            %                      end
            
            %---------------------------更新历史最优y(i)------------------------------
            
            sumEvaluation = sumEvaluation +1;
            
            if fitness(x(MarkOfPar,:),D,roted,MM,train_scale,train_label,test_scale,test_label)>p(MarkOfPar)
                
                p(MarkOfPar)=fitness(x(MarkOfPar,:),D,roted,MM,train_scale,train_label,test_scale,test_label);                     % p(i)为局部最优适应度值；
                
                y(MarkOfPar,:)=x(MarkOfPar,:);                              % y(i)为局部最优位置
                
            end
            %-----------------------------首达次数------------------------------------
            if wheel == 0
                if p(MarkOfPar)<accuracy
                    %                             disp(['第',num2str(runNum),'次运行时间：',num2str(toc)]);   % 程序运行时间
                    successEvaluation = sumEvaluation;
                    accept_iter = interNum;
                    endtime = toc;
                    %                             disp(successEvaluation);
                    %                             disp('success!');
                    wheel = 1;
                    %                             return;
                end
            end
        end
        
        %                         if interNum>100 && mod(interNum,30) ==0
        %                             disp(interNum);
        %                             disp(x);
        %
        %                         end
        %                       if interNum>2000
        %                          disp(['第' num2str(interNum) '代位置矩阵'])
        %                          disp(x)
        %                          disp('=================')
        %                       end
        %
        %                         if mod(interNum,500) == 0
        %                             disp(interNum);
        %                             disp(num);
        %                         end
        
        %==================更新全局最优位置 pg，适应度值 Pbest======================
        for i = 1:N
            if p(i)>fitness(pg,D,roted,MM,train_scale,train_label,test_scale,test_label)
                pg=y(i,:);
            end
            Pbest(interNum)=fitness(pg,D,roted,MM,train_scale,train_label,test_scale,test_label);
        end
        disp([pg '3=====================']);
        %                 disp(['interNum' num2str(interNum) 'Pbest(interNum)' num2str(Pbest(interNum))]);
        
        
        %% ===============================变异===========================-lastIntMinwhere
        sumSwarmFit = 0;
        for parrr = 1:N
            sumSwarmFit = sumSwarmFit + fitness(x(parrr,:),D,roted,MM,train_scale,train_label,test_scale,test_label);
        end
        meanSwarmFit = sumSwarmFit/N;
        
        if interNum>1 && interNum>firstNumMinwhere && interNum<M-lastIntMinwhere && abs(Pbest(interNum)-Pbest(interNum-1))<meanSwarmFit/50
            failSeqNum = failSeqNum+1;
        else
            failSeqNum = 0;
        end
        
        if failSeqNum>=failSeqNumMax
            
            %                      interNum
            num1 = num1+1;
            % 确定变异的粒子和维度
            [MutNum_num,MutDe_num] = MutProMulDem(x,SelMuRat,N,D,NumOfMut,DemOfMut,interNum);
            MutNum_num
            MutDe_num
            %                     [MutNum_num,MutDe_num] = MutPro(x,SelMuRat,N,D,DemOfMut);
            %                    disp(['第' num2str(interNum) '代'])
            %                    disp(['MutNum_num' num2str(MutNum_num)]);
            %                    disp(['MutDe_num' num2str(MutDe_num)]);
            
            for MuNum = 1:NumOfMut
                %% 单个粒子多维同时变异
                xMut=[ x(MutNum_num(MuNum),:); x(MutNum_num(MuNum),:); x(MutNum_num(MuNum),:)];
                for accuNum = 1:3
                    %师兄
                    %                              if rand >0.5
                    %                                 arand = (xMut(accuNum,MutDe_num) +boundary)*((1./(ceil(interNum/1000))).*(1-0.1.^((1-5*(interNum-(floor(interNum/1000))*1000)/5000).^1)));
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) - arand;
                    %                             else
                    %                                 arand =(boundary-xMut(accuNum,MutDe_num))*((1./(ceil(interNum/1000))).*(1-0.1.^((1-5*(interNum-(floor(interNum/1000))*1000)/5000).^1)));
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) + arand;
                    %                              end
                    
                    %在后面人为提高
                    %                               if rand >0.5
                    %                                 arand = (xMut(accuNum,MutDe_num) +boundary)*((1-0.1.^((1-interNum/5000).^1))+0.01*(floor(interNum/500)));
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) - arand;
                    %                             else
                    %                                 arand =(boundary-xMut(accuNum,MutDe_num))*((1-0.1.^((1-interNum/5000).^1))+0.01*(floor(interNum/500)));
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) + arand;
                    %                               end
                    
                    % %                              %反复来回的搜索
                    % %
                    %                               if rand >0.5
                    %                                 arand = (xMut(accuNum,MutDe_num) +boundary)*( 1-0.1.^((1-(interNum-100*floor(interNum/100))/100).^0.8));
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) - arand;
                    %                             else
                    %                                 arand =(boundary-xMut(accuNum,MutDe_num))*( 1-0.1.^((1-(interNum-100*floor(interNum/100))/100).^0.8));
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) + arand;
                    %                              end
                    %                                胡老师策略
                    boundary(MutDe_num,2)
                    xMut(accuNum,MutDe_num(MuNum,:))
                     size(boundary(MutDe_num,2))
                    size(xMut(accuNum,MutDe_num(MuNum,:)))
                    
                    if rand >rand11
                        if rand >0.5
                            if MutDe_num(1) == [1,1]
                            arand = (xMut(accuNum,MutDe_num(MuNum,:)) +cmin )*(1- 0.2.^((1-interNum/M).^1));
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) - arand;
                            MuDistance = [MuDistance arand];
                            end
                            if MutDe_num(1) == [2,1]
                            arand = (xMut(accuNum,MutDe_num(MuNum,:)) +gmin )*(1- 0.2.^((1-interNum/M).^1));
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) - arand;
                            MuDistance = [MuDistance arand];
                            end
                        else
                            if MutDe_num(1) == [1,2]
                            arand =(cmax-xMut(accuNum,MutDe_num(MuNum,:)))*( 1-0.2.^((1-interNum/M).^1));
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) + arand;
                            MuDistance = [MuDistance arand];
                            end
                            if MutDe_num(1) == [2,2]
                            arand =(gmax-xMut(accuNum,MutDe_num(MuNum,:)))*( 1-0.2.^((1-interNum/M).^1));
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) + arand;
                            MuDistance = [MuDistance arand];
                            end
                        end
                    else
                        if rand >0.5
                             if MutDe_num(1) == [1,1]
                             arand = (xMut(accuNum,MutDe_num(MuNum,:)) +cmin)*( 0.2.^((1-interNum/M).^1)-0.1);
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) - arand;
                            MuDistance = [MuDistance arand];
                            end
                            if MutDe_num(1) == [2,1]
                            arand = (xMut(accuNum,MutDe_num(MuNum,:)) +gmin)*( 0.2.^((1-interNum/M).^1)-0.1);
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) - arand;
                            MuDistance = [MuDistance arand];
                            end
                           
                        else
                             if MutDe_num(1) == [1,2]
                           arand =(cmax-xMut(accuNum,MutDe_num(MuNum,:)))*( 0.2.^((1-interNum/M).^1)-0.1);
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) + arand;
                            MuDistance = [MuDistance arand];
                            end
                            if MutDe_num(1) == [2,2]
                           arand =(gmax-xMut(accuNum,MutDe_num(MuNum,:)))*( 0.2.^((1-interNum/M).^1)-0.1);
                            xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) + arand;
                            MuDistance = [MuDistance arand];
                            end
                            
                        end
                    end
                    
                    
                    
                    %% 综合策略
                    % if rand > 0.5
                    
                    %                              if rand >0.5
                    %                                 arand = (xMut(accuNum,MutDe_num(MuNum,:)) +boundary)*((0.000000036802976*interNum.^2 -0.000372022321124*interNum+0.950037201864085).*(sin(pi*(interNum+1290)/70)/2+0.5));
                    %                                 xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) - arand;
                    %                             else
                    %                                 arand =(boundary-xMut(accuNum,MutDe_num(MuNum,:)))*((0.000000036802976*interNum.^2 -0.000372022321124*interNum+0.950037201864085).*(sin(pi*(interNum+1290)/70)/2+0.5));
                    %                                 xMut(accuNum,MutDe_num(MuNum,:)) = xMut(accuNum,MutDe_num(MuNum,:)) + arand;
                    %                              end
                    % else
                    %                              if rand >0.5
                    %                                 arand = (xMut(accuNum,MutDe_num) +boundary)*( 0.1.^((1-interNum/M).^1)-0.1);
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) - arand;
                    %                             else
                    %                                 arand =(boundary-xMut(accuNum,MutDe_num))*( 0.1.^((1-interNum/M).^1)-0.1);
                    %                                 xMut(accuNum,MutDe_num) = xMut(accuNum,MutDe_num) + arand;
                    %                              end
                    % end
                    % for muDem = 1:DemOfMut
                    %                              if rand >0.5
                    %                                 arand = (xMut(accuNum,MutDe_num(muDem)) +boundary)*( 1-0.1.^((1-interNum/M).^1));
                    %                                 xMut(accuNum,MutDe_num(muDem)) = xMut(accuNum,MutDe_num(muDem)) - arand;
                    %                             else
                    %                                 arand =(boundary-xMut(accuNum,MutDe_num(muDem)))*( 1-0.1.^((1-interNum/M).^1));
                    %                                 xMut(accuNum,MutDe_num(muDem)) = xMut(accuNum,MutDe_num(muDem)) + arand;
                    %                             end
                    % end
                    
                    %                             if xMut(accuNum,MutDe_num)>boundary
                    %                                 disp('%%%%%%%%%%%%%%%%%%%%%%%%%')
                    %                                 eee=eee +1;
                    %                             end
                end
                besXmut = xMut(1,:);
                for AccuNum = 2:3
                    if fitness(xMut(AccuNum,:),D,roted,MM,train_scale,train_label,test_scale,test_label)<fitness(besXmut,D,roted,MM,train_scale,train_label,test_scale,test_label)
                        besXmut = xMut(AccuNum,:);
                    end
                end
                %-------------------变异后提高超过 群体平均适应度 了则替换原来的x---------------------
                %                             if fitness(besXmut,D) < fitness(x(MutNum_num(MuNum),:),D)+2
                %                             if fitness(besXmut,D) <meanSwarmFit
                %                                 num11 = num11+1;
                if rand >rand22
                    x(MutNum_num(MuNum),:) = besXmut;
                    sumEvaluation = sumEvaluation +1;
                end
                %                             elseif rand>0.5
                %                                 x(MutNum_num(MuNum),:) = besXmut;
                %                                 num11 = num11+1;
                %                                 反复来回的搜索
                %                             else
                %                               &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
                %% 如果没提高，则向左右跳。但结果不好
                %                                 if rand >0.5
                %                                     x(MutNum_num(MuNum),MutDe_num) = x(MutNum_num(MuNum),MutDe_num)- (x(MutNum_num(MuNum),MutDe_num)-min(x(:,MutDe_num ))) *rand;
                %                                 else
                %                                     x(MutNum_num(MuNum),MutDe_num) =  x(MutNum_num(MuNum),MutDe_num)+ (max(x(:,MutDe_num ))-x(MutNum_num(MuNum),MutDe_num)) *rand;
                %                                 end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % % %                                 如果不好则重复一次不均匀变异，效果比上一个好。
                %                               if rand >0.5
                %                                 arand = ( x(MutNum_num(MuNum),MutDe_num(MuNum,:)) +boundary)*( 0.1.^((1-interNum/M).^1)-0.1);
                %                                  x(MutNum_num(MuNum),MutDe_num(MuNum,:)) =  x(MutNum_num(MuNum),MutDe_num(MuNum,:)) - arand;
                %                             else
                %                                 arand =(boundary- x(MutNum_num(MuNum),MutDe_num(MuNum,:)))*( 0.1.^((1-interNum/M).^1)-0.1);
                %                                  x(MutNum_num(MuNum),MutDe_num(MuNum,:)) =  x(MutNum_num(MuNum),MutDe_num(MuNum,:)) + arand;
                %                              end
                %                                 if fitness(x(MutNum_num(MuNum),:),D)<meanSwarmFit
                %                                     num1111 = num1111+1;
                %                                 end
                %                         elseif rand <0.3
                %                                  arand = rand*(2*boundary-(max(x(:,MutDem_num(MuDemNum)))-min(x(:,MutDem_num(MuDemNum)))));
                %         %                  -------------自己的策略 向外（每一维被选择的概率与距离成反比）-----------
                %                                 if arand >0 && arand<= (min(x(:,MutDem_num(MuDemNum)))+boundary)
                %                                     x(MutNum_num(MuDemNum),MutDem_num(MuDemNum)) = -boundary+arand;
                %                                 else
                %                                     x(MutNum_num(MuDemNum),MutDem_num(MuDemNum)) = -boundary+arand+(max(x(:,MutDem_num(MuDemNum)))-min(x(:,MutDem_num(MuDemNum))));
                %                                 end
                
                %                             end
            end
            
            
            
            failSeqNum = 0;
            MutEither = 1;
        else
            MutEither = 0;
        end
        
        
        %--------------------------前若干代有minwhere---------------------------- interNum>=M-lastIntMinwhere
        if interNum <firstNumMinwhere || interNum+1>=M-lastIntMinwhere|| mod(interNum+1,eachIntMinwhere)==0
            %                      disp('#$%^%&%^&$^*')
            %                      disp(x)
            %%           regress来更新min_where;
            X = zeros(N,D);
            Y = zeros(N,1);
            for k = 1:N;
                X(k,:) = x(k,:);
                Y(k,:)= fitness(x(k,:),D,roted,MM,train_scale,train_label,test_scale,test_label);
            end
           [ min_where]= regress(X,Y,D,regressX,regressXpca,interNum,N);
           
            if min_where(1) == 0
               min_where(1) = rand;
            end
            
            %pcaMatrix= [pcaMatrix pcamatrix];
            %偏最小二乘求minwhere
%             min_where = plsregres(X,Y);
            %                     disp(min_where)
            
            min_wherewhere(interNum,:) = min_where;
            %                     fitness(min_where,D)
            % --------------------更新标准差------------------------------
            for i = 1:D
                estMatr = abs(x(:,i)-min_where(i));
                MinID = 1;
                MinDis = estMatr(1);
                %                       MinID = find(estMatr==min(estMatr));
                % -----------------这一维哪个粒子离他最近----------------------
                for j = 2:N
                    if estMatr(j)<MinDis
                        MinID = j;
                        MinDis = estMatr(MinID);
                    end
                end
                % --这一维上正态分布四个点 期望-min_where与最近点的中点，方差-最近距离
                NormRnd(1:numNormRandMin_where,i) = normrnd...
                    ((min_where(i)+x(MinID,i))/2, MinDis,numNormRandMin_where,1);
            end
            
            if min_where(1) == 0
               min_where(1) = rand;
            end
            
            min_whereFitness(interNum) = fitness(min_where,D,roted,MM,train_scale,train_label,test_scale,test_label);
            
            for  q = 1: numNormRandMin_where
                
                MinNeib = fitness(NormRnd(q,:),D,roted,MM,train_scale,train_label,test_scale,test_label);
                
                if  MinNeib < min_whereFitness(interNum)
                    
                    min_whereFitness(interNum) = MinNeib;
                    
                    min_where =  NormRnd(q,:);
                    
                end
                
            end

            
            %                     x
            %                     fitness(min_where,D)
            %                     sumEvaluation = sumEvaluation +5;
            
        end
        %% 看每一维的聚集程度
%         if interNum<1000
%                 for dd = 1:D
%                     for pp = 1:N
%                         interval = floor((x(pp,dd)+boundary)/(2*boundary/intervalNum));
%                         clusMatrix(dd,min((interval+1),intervalNum)) =  clusMatrix(dd,min((interval+1),intervalNum))+1 ;
%         
%                     end
%                 end
%                 mkdir('mutaResult',[fileName '无变异无全局最优预测早期位置聚集程度分布矩阵正确的ClusMatrix']);
%                 save(['mutaResult/' fileName '无变异无全局最优预测早期位置聚集程度分布矩阵正确的ClusMatrix/第' num2str(interNum) '代分布矩阵ClusMatrix.mat'],'clusMatrix');
%         end
        
    end
    
    %% ===============M代更新完毕 输出
    disp([pg '4=====================']);
  fitness(pg,D,roted,MM,train_scale,train_label,test_scale,test_label)
    disp('success===========');
    %--------------------------程序运行时间------------------------
    if wheel == 0
        endtime = toc;
    end
    %              eee
    %              num
    %              num1
    % num11
    % num1111
   % fitness(pg,D,roted,MM,train_scale,test_scale)
  pg
  disp('5==========')
    accept_iter
    endtime
    result.runtime = [result.runtime endtime];
    result.x = x;
    result.xm = [result.xm  pg'];
    result.fv = [result.fv fitness(pg,D,roted,MM,train_scale,train_label,test_scale,test_label)];
    result.Pbest = [result.Pbest Pbest'];
    result.accept_iter = [result.accept_iter accept_iter];
    result.successEvaluation = [result.successEvaluation successEvaluation];
    result.firstNumPbest = [result.firstNumPbest Pbest(firstNumMinwhere+1)];
    result.vlast =  v;
    result.pcaMatrix = pcaMatrix;
    result.min_whereFitness = [result.min_whereFitness min_whereFitness'];
    %             result.min_wherewhere = min_wherewhere;
    
    for par = 1:N
        fitnesss(par) = fitness(x(par,:),D,roted,MM,train_scale,train_label,test_scale,test_label) ;
        fitnesss(par)
    end
    result.fitness = fitnesss;
end

%-----------Q次取最优那次的Pbest作为figPbest的第一列-----
minIte = 1;
for i = 2:Q
    if result.fv(i)<result.fv(minIte)
        minIte = i;
    end
end
result.figPbest = result.Pbest(:,minIte);
% save(['E:/MATLAB/R2012a/对比粒子群/新comput result/无旋转/' fileName '维度' num2str(D) 'GP-PSOresult.mat'],'result');
%         save(['mutaResult/20140123SelMuRatDminimaFuntion/' fileName ...
%             num2str(w) 'w' num2str(SelMuRat) 'SelMuRat11regress互补策略100代minwhere间隔50代后50代含局部最优变异为均值与当前粒子值变好则替换否则0.5概率变异多维考虑连续失败5维度' ...
%             num2str(D)  num2str(NumOfMut) '个粒子' num2str(DemOfMut) '维变异c=0.9c3=' num2str(c3) '.mat'],'result');
%        num2str(numOfMut)
% % %% b变异距离统计操作
% MuDistance0 = [];
% MuDistance1 = [];
% MuDistance2 = [];
% for z = 1:length(MuDistance);
%    if mod(z,30)==1 || mod(z,30)==15
%        MuDistance0 = [MuDistance0 MuDistance(z)];
%    end
% end
% for zz = 1:length(MuDistance0);
%    if mod(zz,6)==1 || mod(zz,6)==2
%        MuDistance1 = [MuDistance1 MuDistance0(zz)];
%    end
% end
% for zzz = 1:length(MuDistance1)
%    if mod(zzz,30)==1 || mod(zzz,30)==2
%        MuDistance2 = [MuDistance2 MuDistance1(zzz)];
%    end
% end

% 变异30个粒子30维统计变异距离操作
% MuDistance0 = [];
% for z = 1:length(MuDistance)
%     if mod(z,DemOfMut*NumOfMut*3)<DemOfMut*NumOfMut && mod(z,DemOfMut)==1
%         MuDistance0 = [MuDistance0 MuDistance(z)];
%     end
% end
% 
% save(['mutaResult/其他实验结果/' fileName 'rand22为' num2str(rand22) 'rand11为' num2str(rand11) 'rou为0.2的扩展的不均匀变异适当数量点GP-PSOresult.mat'],'result');
% save(['mutaResult/其他实验结果/' fileName 'rand22为' num2str(rand22) 'rand11为' num2str(rand11) 'rou为0.2的扩展的不均匀变异适当数量点MuDistance2000代rand值.mat'],'MuDistance0');
% save(['mutaResult/其他实验结果/' fileName 'rand22为' num2str(rand22) 'rou为0.2的标准不均匀变异GP-PSOresult.mat'],'result');
% save(['mutaResult/其他实验结果/' fileName 'rand22为' num2str(rand22) 'rou为0.2的标准不均匀变异MuDistance2000代rand值.mat'],'MuDistance0');
% save(['mutaResult/其他实验结果/' fileName num2str(rand11) 'rand11会议论文GP-PSOresult.mat'],'result');
% save(['mutaResult/其他实验结果/' fileName num2str(rand11) 'rand11会议论文MuDistance2000代rand值.mat'],'MuDistance2');
% end

% save('mutaResult/AckleyFunction分布矩阵Aresult/GP-PSOresult.mat','result');

% save(['mutaResult/简单函数无变异pcamatrix/' fileName 'GP-PSOresult.mat'],'result');

% save(['mutaResult/' fileName '无变异无全局最优预测早期位置聚集程度分布矩阵正确的ClusMatrix/GP-PSOresult.mat'],'result');
% save(['mutaResult/测试非均匀变异及扩展后的0.1换为rand的结果/' fileName '标准的不均匀变异randGP-PSOresult.mat'],'result');
% save(['E:/MATLAB/R2012a/对比粒子群/新comput result/加旋转/' fileName '维度' num2str(D) '/GP-PSOresult后跑.mat'],'result');




