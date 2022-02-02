function [beta, betaCurv, betaBest, lambda, lambdaCurv, lambdaBest, ...
    lambdaGridMinNorm] = minNorm_lcurve_bestRegul( G , b, source)

[u,s,v] = csvd(G);

% Find corner (lambda) of the L-curve 
% curvature is computed with lcfun (lambda) or the curvature function (lambdaCurv)
[lambda,rho,eta,lambdaGridMinNorm,lambdaCurv] = l_curve_time(u,s,b,'Tikh');

% Compute the minimum-norm estimates (betas) for all lambdas
betaMinNorm = cell(1, numel(lambdaGridMinNorm));
for i = 1:numel(lambdaGridMinNorm)
    betaMinNorm{i} = zeros(size(G, 2), size(b, 2));
    for ndx = 1:size(b, 2)
        betaMinNorm{i}(: ,ndx) = tikhonov( u , s , v , b(:, ndx), lambdaGridMinNorm(i));
    end
end

% compute beta for lambda
beta = zeros(size(betaMinNorm{1}));
for i = 1:size(b, 2)
    beta(:, i) = tikhonov(u, s, v, b(:, i), lambdaCurv);
end

% compute beta for lambdaCurv
betaCurv = zeros(size(betaMinNorm{1}));
for i = 1:size(b, 2)
    betaCurv(:, i) = tikhonov(u, s, v, b(:, i), lambdaCurv);
end

% compute beta for lambda chosen based on min MSE
% first compute MSE
winERP = 46:180;
mseTruth = zeros(length(betaMinNorm),1);
for iLambda = 1:length(betaMinNorm)
    mseTruth(iLambda) = sum(sum((betaMinNorm{iLambda}-source).^2));
%     [aucAve(iLambda), energyAve(iLambda), mseAve(iLambda)] = computeMetrics(betaMinNorm{iLambda}(:,winERP),source(:,winERP));
end
% then get min MSE to pick the best lambda
lambdaBest = lambdaGridMinNorm(mseTruth==min(mseTruth));
% compute beta with that lambda
betaBest = zeros(size(betaMinNorm{1}));
for i = 1:size(b, 2)
    betaBest(:, i) = tikhonov(u, s, v, b(:, i), lambdaBest);
end

% figure;
% semilogy(lambdaGridMinNorm,mseTruth,'b','linewidth',2)
% hold on
% semilogy([lambda lambda],ylim,'g:','linewidth',2)
% semilogy([lambdaCurv lambdaCurv],ylim,'b:','linewidth',2)
% semilogy([lambdaBest lambdaBest],ylim,'r:','linewidth',2)
% xlabel('lambda');ylabel('log MSE')

% figure;
% semilogy(lambdaGridMinNorm,aucAve,'b','linewidth',2)
% hold on
% semilogy([lambda lambda],ylim,'g:','linewidth',2)
% semilogy([lambdaCurv lambdaCurv],ylim,'b:','linewidth',2)
% semilogy([lambdaBest lambdaBest],ylim,'r:','linewidth',2)
% xlabel('lambda');ylabel('log MSE')

