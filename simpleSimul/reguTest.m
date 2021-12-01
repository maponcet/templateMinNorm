

addpath([pwd filesep 'reguTime' filesep]);

load('matrixMinNorm.mat')

groundTruth = zeros(18,90);
groundTruth(:,46:90) = repmat( [1,1,zeros(1,16)]',1,45);
% groundTruth= repmat( [1,1,zeros(1,16)]',1,45);

noise = 0*randn(size(data));
datan= data+noise;

[betaMinNormBest, betaMinNormOver, bestLambdaOver, lambdaGridOver] = minNorm_lcurve( templateOverlap , datan );
[betaMinNormBestNoOver, betaMinNormNoOver, bestLambdaNoOver, lambdaGridNoOver] = minNorm_lcurve( templateNoOverlap , datan );

[betaMinNormBestGCV, betaMinNormOverGCV, lambdaOver, gcvErrorMinNormOver, lambdaGridOverGCV] = minimum_norm( templateOverlap , datan,10 );
[betaMinNormBestNoOverGCV, betaMinNormNoOverGCV, lambdaNoOver, gcvErrorMinNormNoOver, lambdaGridNoOverGCV] = minimum_norm( templateNoOverlap , datan,10 );

for iLambda = 1:length(betaMinNormOver)
    mseOverlap(iLambda) = sum(sum((betaMinNormOver{iLambda}-groundTruth).^2));
    mseNoOverlap(iLambda) = sum(sum((betaMinNormNoOver{iLambda}-groundTruth).^2));
end


figure;
semilogy(lambdaGridOver,mseOverlap)
hold on
semilogy(lambdaGridNoOver,mseNoOverlap)

semilogy([bestLambdaOver bestLambdaOver],ylim,'k--','linewidth',2)

semilogy([bestLambdaOver bestLambdaOver],ylim,'k:','linewidth',2)

semilogy([bestLambdaOverGCV bestLambdaOverGCV],ylim,'r--','linewidth',2)
semilogy([bestLambdaNoOverGCV bestLambdaNoOverGCV],ylim,'r:','linewidth',2)

ylabel('MSE to ground truth')
xlabel('Lambda')
legend('Overlap','No Overlap','Chosen lambda Overlap','Chosen Lambda No Over','GCV Overlap','GCV No Over')


figure;
subplot(1,2,1);
semilogy(lambdaGridOver,mseOverlap)
hold on;
semilogy(lambdaGridNoOver,mseNoOverlap)
semilogy([bestLambdaOver bestLambdaOver],ylim,'k--','linewidth',2)
semilogy([bestLambdaOver bestLambdaOver],ylim,'k:','linewidth',2)
semilogy([bestLambdaOverGCV bestLambdaOverGCV],ylim,'r--','linewidth',2)
semilogy([bestLambdaNoOverGCV bestLambdaNoOverGCV],ylim,'r:','linewidth',2)
ylabel('MSE to ground truth')
xlabel('Lambda')
subplot(1,2,2); hold on
plot(lambdaGridOverGCV,gcvErrorMinNormOver)
plot(lambdaGridNoOverGCV,gcvErrorMinNormNoOver)
plot([bestLambdaOver bestLambdaOver],ylim,'k--','linewidth',2)
plot([bestLambdaOver bestLambdaOver],ylim,'k:','linewidth',2)
plot([bestLambdaOverGCV bestLambdaOverGCV],ylim,'r--','linewidth',2)
plot([bestLambdaNoOverGCV bestLambdaNoOverGCV],ylim,'r:','linewidth',2)
ylabel('gcv error')
xlabel('Lambda')
legend('Overlap','No Overlap','Chosen lambda Overlap','Chosen Lambda No Over','GCV Overlap','GCV No Over')


tic
[betaMinNormBest, betaMinNormOver, bestLambdaOver, lambdaGridOver] = minNorm_lcurve( templateOverlap , datan );
toc
tic
[betaMinNormBest, lambdaBest] = minNormFast_lcurve( templateOverlap , datan );
toc
