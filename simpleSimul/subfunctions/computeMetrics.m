function [auc, energy,mse] = computeMetrics(beta,indSources,sourceOverTime)
% beta is the retrieved signal (beta values output from minimum_norm) per
% ROI over time (here 18 ROI x 90 timepoints)
% indSources is a vector of the activated source indexes (so length should
% be equal to the total number of activated sources)
% sourceOverTime is the pure (no noise) simulated activity for each ROI
% over time

%%%%% AUC
tmp = zeros(size(beta,1),1);
tmp(indSources,:) = 1;
aucTime = zeros(1,size(beta,2));
for nT = 1:size(beta,2)
    aucTime(nT) = rocArea( abs(beta(:,nT)) , tmp );
end
auc = mean(aucTime);

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
    % the 2 lines below are equal
%     mseBeta(nT) = sum( (normSource(indSources) - norm_beta(indSources,nT)).^2 ) / numel(normSource(indSources));
    mseBeta(nT) = norm(normSource(indSources) - norm_beta(indSources,nT),'fro').^2 / numel(norm_beta(indSources,nT));
end
mse = mean(mseBeta);

% % test frobenius:
% a = randn(10,1); b=randn(10,1);
% norm(a-b,'fro').^2
% sum((a-b).^2)