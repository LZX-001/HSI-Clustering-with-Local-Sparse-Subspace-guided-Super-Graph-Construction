%--------------------------------------------------------------------------
% This function takes a DxN matrix of N data points in a D-dimensional 
% space and returns a NxN coefficient matrix of the sparse representation 
% of each data point in terms of the rest of the points
% Y: DxN data matrix
% affine: true if enforcing the affine constraint, false otherwise
% thr1: stopping threshold for the coefficient error ||Z-C||
% thr2: stopping threshold for the linear system error ||Y-YZ||
% maxIter: maximum number of iterations of ADMM
% C2: NxN sparse coefficient matrix
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function C2 = admmOutlier_mat_func(Y,Con, affine, alpha, thr, maxIter)

if (nargin < 3)
    % default subspaces are linear
    affine = false; 
end
if (nargin < 4)
    % default regularizarion parameters
    alpha = 20;
end
if (nargin < 5)
    % default coefficient error threshold to stop ALM
    % default linear system error threshold to stop ALM
    thr = 2*10^-4; 
end
if (nargin < 6)
    % default maximum number of iterations of ALM
    maxIter = 200; 
end

if (length(alpha) == 1)
    alpha1 = alpha(1);
    alpha2 = alpha(1);
    alpha3 = alpha(1);
elseif (length(alpha) == 2)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    alpha3 = alpha(2);
elseif (length(alpha) == 3)
    alpha1 = alpha(1);
    alpha2 = alpha(2);
    alpha3 = alpha(3);
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

if (~affine)
    % initialization
    A = inv(mu1 * (P' * P) + mu2 * eye(N + D));
    C1 = zeros(N + D, N);
    Lambda1 = zeros(D, N);
    Lambda2 = zeros(N + D, N);
    err1 = 10 * thr1; err2 = 10 * thr2;
    ii = 1;
    % ADMM iterations
    while ((err1(ii) > thr1 ||err2(ii)>thr2) && ii < maxIter)
        % updating Z
        Z = A * (mu1 * P' * (Y + Lambda1 / mu1) + mu2 * (C1 - Lambda2 / mu2));
        Z(1:N, :) = Z(1:N, :) - diag(diag(Z(1:N, :)));
        % updating C
        C2 = max(0, (abs(Z + Lambda2 / mu2) - 1 / mu2 * ones(N + D, N))) .* sign(Z + Lambda2 / mu2);
        C2(1:N, :) = C2(1:N, :) - diag(diag(C2(1:N, :)));
        C2(1:N,:)=  C2(1:N,:).*Con;%主要的，做邻居关系限制的
        % updating Lagrange multipliers
        Lambda1 = Lambda1 + mu1 * (Y - P * Z);
        Lambda2 = Lambda2 + mu2 * (Z - C2);
        % computing errors
        err1(ii + 1) = errorCoef(Z, C2);
        err2(ii + 1) = errorLinSys(P, Z);
        disp([err1(ii+1) err2(ii+1)]);
        %
        C1 = C2;
        ii = ii + 1;
    end
    % fprintf('err1: %2.4f, err2: %2.4f, iter: %3.0f \n', err1(end), err2(end), ii);
else
    % initialization
    delta = [ones(N, 1); zeros(D, 1)];
    A = inv(mu1 * (P' * P) + mu2 * eye(N + D) + mu2 * (delta * delta'));
    C1 = zeros(N + D, N);
    Lambda1 = zeros(D, N);
    Lambda2 = zeros(N + D, N);
    lambda3 = zeros(1, N);
    err1 = 10 * thr1; err2 = 10 * thr2; err3 = 10 * thr1;
    ii = 1;
    % ADMM iterations
    while ((err1(ii) > thr1 || err2(ii) > thr2 || err3(ii) > thr1) && ii < maxIter)
        % updating Z
        Z = A * (mu1 * P' * (Y + Lambda1 / mu1) + mu2 * (C1 - Lambda2 / mu2) + mu2 * delta * (ones(1, N) - lambda3 / mu2));
        Z(1:N, :) = Z(1:N, :) - diag(diag(Z(1:N, :)));
        % updating C
        C2 = max(0, (abs(Z + Lambda2 / mu2) - 1 / mu2 * ones(N + D, N))) .* sign( Z + Lambda2 / mu2);
        C2(1:N, :) = C2(1:N, :) - diag(diag(C2(1:N, :)));
        % updating Lagrange multipliers
        Lambda1 = Lambda1 + mu1 * (Y - P * Z);
        Lambda2 = Lambda2 + mu2 * (Z - C2);
        lambda3 = lambda3 + mu2 * (delta'*Z - ones(1, N));
        % computing errors
        err1(ii + 1) = errorCoef(Z, C2);
        err2(ii + 1) = errorLinSys(P, Z);
        err3(ii + 1) = errorCoef(delta' * Z, ones(1, N));
        %
        C1 = C2;
        ii = ii + 1;
    end
    % fprintf('err1: %2.4f, err2: %2.4f, err3: %2.4f, iter: %3.0f \n', err1(end), err2(end), err3(end), ii);
end