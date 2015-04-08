function x = LHspso(N,D)
% boundary = 5.12;
% N = 40;
% D = 30;

D=2;
cmin=0.001;
cmax=1;
gmin=0.001;
gmax=10;

%--------区间长度------
%manus = (2*boundary)/N;
manus = [(cmax-cmin)/N;(gmax-gmin)/N];
%--------生成区间------
internle = [];
matrix = [];
x = zeros(N,D);
%{

for i = 1:N
    internle(i,1) = -boundary+(i-1)*manus;
    internle(i,2) = -boundary+i*manus;
end
%}
% disp(internle);
% size(internle);
%--------每段区间随机取30个数,生成矩阵matrix-----
for i = 1:N
   % matrix(:,i) = rand(1,1)*manus-[(cmax-cmin);(gmax-gmin)]+(i-1)*manus;
    matrix(:,i) = [(cmax-cmin);(gmax-gmin)]*rand+manus;
end

for m=1:N
    if matrix(1,m)>cmax
        matrix(1,m) = cmax*rand;
    end
    if matrix(1,m)<cmin
        matrix(1,m) = cmin+rand;
    end
    
    if matrix(2,m)>gmax
        matrix(2,m) = gmax*rand;
    end
    if matrix(2,m)<gmin
        matrix(2,m) = gmin+rand;
    end
end

% disp(matrix);
%-------生成随机点-----------------
for j = 1:N
    
    for i = 1:D
        siz = size(matrix);
        
        idx = randperm(siz(2));
        matrix(i,:) = matrix(i,idx);
    end
    
     x (j,:) = matrix(:,1);
    matrix = matrix(:,2:siz(2));
  
end


