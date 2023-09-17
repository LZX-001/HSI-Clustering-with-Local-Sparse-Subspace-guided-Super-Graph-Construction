function C2 = graphlearning(Y,Con, alpha1,alpha2, thr, maxIter)
if (nargin < 5)
    % default coefficient error threshold to stop ALM
    % default linear system error threshold to stop ALM
    thr = 2*10^-4; 
end
if (nargin < 6)
    % default maximum number of iterations of ALM
    maxIter = 200; 
end

if (length(thr) == 1)
    thr1 = thr(1);
    thr2 = thr(1);
elseif (length(thr) == 2)
    thr1 = thr(1);
    thr2 = thr(2);
end

[D, N] = size(Y);

gamma = alpha3 / norm(Y, 1);
P = [Y, eye(D)/gamma];

% setting penalty parameters for the ADMM
mu1 = alpha1 * 1 / computeLambda_mat(Y, P);
mu2 = alpha2 * 1;


% initialization
A = inv(mu1 * (P' * P) + mu2 * eye(N + D));
C1 = zeros(N + D, N);
Lambda1 = zeros(D, N);
Lambda2 = zeros(N + D, N);
err1 = 10 * thr1; err2 = 10 * thr2;
ii = 1;
% ADMM iterations
while ((err1(ii) > thr1 || err2(ii) > thr2) && ii < maxIter)
    %updating Z
    Z = A * (mu1 * P' * (Y + Lambda1 / mu1) + mu2 * (C1 - Lambda2 / mu2));
    Z(1:N, :) = Z(1:N, :) - diag(diag(Z(1:N, :)));
    updating C
    C2 = max(0, (abs(Z + Lambda2 / mu2) - 1 / mu2 * ones(N + D, N))) .* sign(Z + Lambda2 / mu2);
    C2(1:N, :) = C2(1:N, :) - diag(diag(C2(1:N, :)));
    
    %updating Lagrange multipliers
    Lambda1 = Lambda1 + mu1 * (Y - P * Z);
    Lambda2 = Lambda2 + mu2 * (Z - C2);
    %computing errors
    err1(ii + 1) = errorCoef(Z, C2);
    err2(ii + 1) = errorLinSys(P, Z);
    disp([err1(ii+1) err2(ii+1)]);
    C1 = C2;
    ii = ii + 1;
end


