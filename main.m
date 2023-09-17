clear;
addpath('.\SSC_ADMM_v1.1');
addpath('.\ClusteringMeasure');
addpath('.\EntropyRateSuperpixel-master');
Dataset='PU'%'IP','SV','PU'
if(strcmp(Dataset,'IP'))
    R=importdata('.\data\Indian_pines_corrected.mat');
    gt=importdata('.\data\Indian_pines_gt.mat');
    [m,n,d]=size(R);
    rho=0.3;%0.2
    lambda1=25;
    lambda2=0.1;
    knearest=10;
    r=60;
    K=75;
end
if(strcmp(Dataset,'SA'))
    R=importdata('.\data\Salinas_corrected.mat');
    gt=importdata('.\data\Salinas_gt.mat');
    [m,n,d]=size(R);
    rho=0.2;%0.2
    lambda1=250;
    lambda2=10;
    knearest=10;
    r=50;
    K=70;
end
if(strcmp(Dataset,'PU'))
    R=importdata('.\data\PaviaU_corrected.mat');
    gt=importdata('.\data\PaviaU_gt.mat');
    [m,n,d]=size(R);
    rho=0.3;%0.2
    lambda1=15;
    lambda2=10;
    knearest=10;
    r=60;
    K=70;
end
tic
X1=reshape(R,m*n,d);

label=generateSuperpixel(double(R),K);
for i=1:K
    f=find(label==i);
    Y(i,:)=mean(X1(f,:));
end

M1=createGraph(label,m,n);%neighbor
    
num_class=max(gt(:));
[ssc_coef,A,distance]=self_method(Y,M1,lambda1,lambda2,knearest,r,rho);
gt1=reshape(gt,m*n,1);
f=find(gt1==0);
gt1(f)=[];
re=zeros(m*n,1);
for i=1:10
    result = spectral_clustering(A, double(num_class));
    for j=1:K
            f1=find(label==j);
            re(f1)=result(j);
     end
     rr(i,:)=re;
     re(f)=[];
     [acc1(i) acc2(i) kappa(i) F1(i) P(i) Recall(i) Purity(i) ARI(i) DX(i,:) c(i,:)]=Accuracy(re,gt1);
     
    [AX nmi(i) avgent] = compute_nmi(gt1,re);
end
[max_acc1 index]=max(acc1);
D_max=DX(index,:);
c_max=c(index,:);
%r_max=rr(index,:);
%save('.\S2_max.mat','D_max');
disp('avg acc, std and max acc');
disp([mean(acc1) std(acc1) max_acc1]);
disp('avg AA, std and max AA');
disp([mean(acc2) std(acc2) acc2(index)]);
disp('avg kappa, std and max kappa');
disp([mean(kappa) std(kappa) kappa(index)]);
disp('avg nmi, std and max nmi');
disp([mean(nmi) std(nmi) nmi(index)]);
disp('avg F1, std and max F1');
disp([mean(F1) std(F1) F1(index)]);
disp('avg P, std and max P');
disp([mean(P) std(P) P(index)]);
disp('avg Recall, std and max Recall');
disp([mean(Recall) std(Recall) Recall(index)]);
disp('avg Purity, std and max Purity');
disp([mean(Purity) std(Purity) Purity(index)]);
disp('avg ARI, std and max ARI');
disp([mean(ARI) std(ARI) ARI(index)]);

toc


