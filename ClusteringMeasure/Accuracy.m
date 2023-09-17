function [OA,AA,KAPPA,F1,P,R,Purity,ARI,C,c] = Accuracy(C,gt)
 [C,c] = bestMap(gt,C);
[F1,P,R]=compute_f(gt,C);
Purity=purity(C,gt);
ARI=RandIndex(C,gt);

 OA = length(find(gt == C))/length(gt);
len=max(gt(:));
x1=[];
x2=[];
x3=[];
for ii=1:len
    x1(ii)=length(find(gt==ii));
    x2(ii)=length(find(C==ii));
    x3(ii)=length(find(gt == C&gt==ii))/x1(ii);
end
%disp(x3);
AA=mean(x3);
pe=sum(x1.*x2)/(length(gt)^2);
KAPPA=(OA-pe)/(1-pe);
end