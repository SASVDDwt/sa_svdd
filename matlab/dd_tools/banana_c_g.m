function [  ] = banana_c_g( a,n,c,g,pic_num )

%{
a=train_scale;
a = oc_set(a,'1');
a = target_class(a);
figure(pic_num); hold on;
h = scatterd(a);

fracrej = 1/(n*c);
w2 = svdd(a,fracrej,g);
h = plotc(w2,'r--');hold on;
option= ['-s 5',' -g ',num2str(g),' -c ',num2str(c)];
model=svmtrain(train_label,train_scale,option);
plot(model.SVs(:,1),model.SVs(:,2),'ko');hold on;
%}

%a = gendatb([100 0]);
%a = oc_set(a,'1');
%a = target_class(a);

figure(pic_num);
hold on;

h = scatterd(a);
fracrej = 1/(n*c);
w2 = svdd(a,fracrej,g);
h = plotc(w2,'r--');

end

