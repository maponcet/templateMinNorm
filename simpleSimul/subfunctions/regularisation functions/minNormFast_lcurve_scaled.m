function [betaMinNormBest, lambda] = minNormFast_lcurve_scaled( G , b)

[u,s,v] = csvd(G);
lambda = l_curve_time(u,s,b,'Tikh');

betaMinNormBest = zeros([size(G,2) size(b, 2)]);
for ll = 1:size(b, 2)
    betaMinNormBest(:, ll) = tikhonov(u, s, v, b(:, ll), lambda);
end

fitData = G*betaMinNormBest;
scalingCoef = fitData(:)\b(:);
betaMinNormBest = betaMinNormBest*scalingCoef;
