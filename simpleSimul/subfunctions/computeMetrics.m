function [auc, energy,mse] = computeMetrics(beta,indSources,sourceOverTime)
% compute mean over time for AUC, relative energy, MSE
% the input is the average (all sbj) retrieved signal per ROI over time =
% 18 ROI * 90 timepoints matrix of the output betas from minimum_norm
% function

%%%%% AUC
tmp = zeros(size(beta,1),1);
tmp(indSources,:) = 1;
% instead of using 1 and 0 in the matrix, use the real amount of simulated
% activity? Nooo rocArea only uses 1 and 0s
% AUC computed without normalising: it wouldn't change the results
aucTime = zeros(1,size(beta,2));
for nT = 1:size(beta,2)
    aucTime(nT) = rocArea( abs(beta(:,nT)) , tmp );
end
auc = mean(aucTime); % auc calculated from the average activity (beta)

%%%%% relative energy
norm_beta= zeros(size(beta,1),size(beta,2)); relEnergy= zeros(1,size(beta,2));
for nT = 1:size(beta,2)
    norm_beta(:,nT) = beta(:,nT) / max( abs(beta(:,nT)) ); % normalise estimated sources
    relEnergy(nT) = sum( abs( norm_beta(indSources,nT) ) ) / sum( abs(norm_beta(:,nT)) );
end
energy = mean(relEnergy);

%%%%% mse
% paper: norm(sourceValOverTime(:,nT)-betaReg(:,nT),'fro').^2 / numel(betaReg(:,nT))
mseBeta= zeros(1,size(beta,2));
for nT = 1:size(beta,2)
    normSource = sourceOverTime(:,nT) / max(abs(sourceOverTime(:,nT))); % normalise simulated source
    mseBeta(nT) = sum( (normSource(indSources) - norm_beta(indSources,nT)).^2 ) / sum( (normSource(indSources)).^2 );
end
mse = mean(mseBeta);

