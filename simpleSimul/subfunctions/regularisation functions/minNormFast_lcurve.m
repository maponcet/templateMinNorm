function [betaMinNormBest, lambda,residualNorm,solutionNorm,regulOK] = minNormFast_lcurve( G , b)

[u,s,v] = csvd(G);
[lambda,residualNorm,solutionNorm,reg_param] = l_curve_modified(u,s,b,'Tikh');

betaMinNormBest = zeros([size(G,2) size(b, 2)]);
for ll = 1:size(b, 2)
    betaMinNormBest(:, ll) = tikhonov(u, s, v, b(:, ll), lambda);
end


if solutionNorm > 300 || sum(lambda>reg_param) < 3 || sum(lambda<reg_param)<3
    regulOK = 0;
else
    regulOK = 1;
end