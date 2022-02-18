function [auc, energy,mseNorm] = computeMetrics(beta,sourceOverTime)
% beta is the retrieved signal (beta values output from minimum_norm) per
% ROI over time (here 18 ROI x 90 timepoints)
% sourceOverTime is the pure (no noise) simulated activity for each ROI
% over time

%%%%% AUC
aucTime = zeros(1,size(beta,2));
srcOn = sourceOverTime~=0;
if isempty(find(srcOn~=1))
    auc=NaN;
else
    for nT = 1:size(beta,2)
        aucTime(nT) = rocArea( abs(beta(:,nT)) , srcOn(:,nT) );
    end
    auc = mean(aucTime);
end

%%%%% relative energy
norm_beta= zeros(size(beta,1),size(beta,2)); relEnergy= zeros(1,size(beta,2));
for nT = 1:size(beta,2)
    activ_sources = sourceOverTime(:,nT)~=0;
    norm_beta(:,nT) = beta(:,nT) / max( abs(beta(:,nT)) ); % normalise estimated sources
    relEnergy(nT) = sum( abs( norm_beta(activ_sources,nT) ) ) / sum( abs(norm_beta(:,nT)) );
end
energy = mean(relEnergy);

%%%%% mse
% paper: norm(sourceValOverTime(:,nT)-betaReg(:,nT),'fro').^2 / numel(betaReg(:,nT))
mseBeta= zeros(1,size(beta,2));
mseBetaNorm= zeros(1,size(beta,2));
for nT = 1:size(beta,2)
    normSource = sourceOverTime(:,nT) / max(abs(sourceOverTime(:,nT))); % normalise simulated source
    % the 2 lines below are equal
%     mseBeta(nT) = sum( (normSource - norm_beta(:,nT)).^2 ) / numel(normSource);
    mseBetaNorm(nT) = norm(normSource - norm_beta(:,nT),'fro').^2 / numel(norm_beta(:,nT));
    mseBeta(nT) = norm(sourceOverTime(:,nT) - beta(:,nT),'fro').^2 / numel(beta(:,nT));
end
mseNorm = mean(mseBetaNorm); 
mse = mean(mseBeta); % not meaningful

% % test frobenius:
% a = randn(10,1); b=randn(10,1);
% norm(a-b,'fro').^2
% sum((a-b).^2)