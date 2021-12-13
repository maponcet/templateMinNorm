function [betaMinNormBest, lambda] = minNormFast( G , b, nLambda )

% Introduction of the prior on the sources (if no prior recquired, R_sq is
% the identity).
%G = G * R_sq;
% Optimization of the regularization parameter according to the gcv error
[u,s,v] = csvd(G);
[lambda] = gcv(u,s,b,'Tikh', nLambda);

betaMinNormBest = zeros([size(G,2) size(b, 2)]);
for i = 1:size(b, 2)
    betaMinNormBest(:, i) = tikhonov(u, s, v, b(:, i), lambda);
end

fitData = G*betaMinNormBest;

scalingCoef = fitData(:)\b(:);
betaMinNormBest = betaMinNormBest*scalingCoef;
