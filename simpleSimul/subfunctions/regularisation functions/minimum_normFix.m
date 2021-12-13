function [betaMinNormBest] = minimum_normFix( model , data, nLambda )
%[betaMinNormBest, rho, eta] = minimum_normFix( model , data, nLambda )

% Introduction of the prior on the sources (if no prior recquired, R_sq is
% the identity).
%G = G * R_sq;
% Optimization of the regularization parameter according to the gcv error
[u,s,v] = csvd(model);
rho = zeros(size(data,2),1);
eta = zeros(size(data,2),1);
betaMinNormBest = zeros(size(model,2),size(data,2));
for tt = 1:size(data, 2)
    [betaMinNormBest(:, tt),rho(tt),eta(tt)] = tikhonov(u, s, v, data(:, tt), nLambda);
end
