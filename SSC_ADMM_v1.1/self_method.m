function [C1,A,distance] = self_method(Y,Con,alpha1,alpha2,k,r,rho)
Y=DataProjection(Y',r);
%Y=(Y-min(min(Y)))/(max(max(Y))-min(min(Y)));
[D, N] = size(Y);
maxIter = 200; %
thr2=2*10^-4;
thr1=thr2; 
% setting penalty parameters for the ADMM
mu1 = alpha1 * 1 / computeLambda_mat(Y, Y);
mu2 = alpha1 * 1;
mu3=alpha1;
mu4=alpha2;
T = inv(mu1 * (Y' * Y) + mu3 * eye(N));%+mu4*eye(N));
C1 = zeros(N, N);
E= zeros(D, N);
A=zeros(N,N);
Lambda1 = zeros(N, N);
Lambda2=zeros(N,N);
err1 = 10 * thr1; err2 = 10 * thr2;
ii = 1;
% ADMM iterations
while ((err1(ii) > thr1||err2(ii)>thr2 ) && ii < maxIter)

    % updating Z
    Z = T * (mu1 * Y' * (Y - E) + mu3 * (C1 - Lambda1 / mu3));%+mu4*((Y'*W)*(Y'*W)'));
    Z = Z - diag(diag(Z));
    % updating C
    C2 = max(0, (abs(Z + Lambda1 / mu3) - 1 / mu3 * ones(N , N))) .* sign(Z + Lambda1 / mu3);  
    C2 = C2 - diag(diag(C2));

    % updating E
    E = max(0, (abs(Y-Y*Z) - mu2 / mu1 * ones(D, N))) .* sign(Y-Y*Z);
    % constructing L
    Dt=sum(A,1);
    Dm=diag(Dt);
    L=Dm-A;
    %updating X
    X=(mu3*C2-Lambda2)*inv(2*mu4*L+mu3*eye(N,N));
    X = X - diag(diag(X));
    %updating A
    distance=(pdist2(X',X').^2)/4;
    for i=1:N
        temp=sort(distance(i,:));
        A(:,i)=max(0,(1+sum(temp(1:k)))/k-distance(:,i));
    end
    A=(A+A')/2;
    A=A.*Con;
    Lambda1 = Lambda1 + mu3 * (Z - C2);
    Lambda2=Lambda2+mu3*(X-C2);
    % computing errors
    err1(ii + 1) = errorCoef(Z, C2);
    err2(ii + 1) = errorCoef(X,C2);
    disp([err1(ii+1) err2(ii+1)]);
    %disp(errorCoef((Y'*W)*(Y'*W)',C2))
    %
    C1 = C2;
    %C1 = BuildAdjacency(thrC(C1, rho));
    
    ii = ii + 1;
end
C1 = BuildAdjacency(thrC(C1, rho));
A=BuildAdjacency(thrC(A, rho));
end
