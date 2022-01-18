function [betaMinNormBest, betaMinNorm, lambda, lambdaGridMinNorm,lambdaCurv] = minNorm_lcurve( G , b)

% Introduction of the prior on the sources (if no prior recquired, R_sq is
% the identity).
%G = G * R_sq;
% Optimization of the regularization parameter according to the gcv error
[u,s,v] = csvd(G);
%[lambda, gcvErrorMinNorm, lambdaGridMinNorm] = gcv(u,s,b,'Tikh', nLambda);

[lambda,rho,eta,lambdaGridMinNorm,lambdaCurv] = l_curve_time(u,s,b,'Tikh');

% Compute the minimum-norm estimates based on the gcv error.
betaMinNorm = cell(1, numel(lambdaGridMinNorm));
for i = 1:numel(lambdaGridMinNorm)
    betaMinNorm{i} = zeros(size(G, 2), size(b, 2));
    for ndx = 1:size(b, 2)
        betaMinNorm{i}(: ,ndx) = tikhonov( u , s , v , b(:, ndx), lambdaGridMinNorm(i));
    end
end

betaMinNormBest = zeros(size(betaMinNorm{1}));
for i = 1:size(b, 2)
    betaMinNormBest(:, i) = tikhonov(u, s, v, b(:, i), lambdaCurv);
end


