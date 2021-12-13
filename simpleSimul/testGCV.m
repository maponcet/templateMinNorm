%%%%%%%%%%%%%%%%%%%%
%%%% test gcv, try to understand why gcv does not go up
% plot data, fit, residual for different lambda
% minimum_norm(template, data, 10) = minimum_norm( G , b, nLambda )

addpath('/Users/marleneponcet/Documents/Git/svndl_code/alesToolbox');
addpath([pwd filesep 'subfunctions' filesep]);

load matrixMinNorm.mat

nLambda = 20;

% Optimization of the regularization parameter according to the gcv error
[u,s,v] = csvd(templateNoOverlap);
[uO,sO,vO] = csvd(templateOverlap);
[lambda, gcvErrorMinNorm, lambdaGridMinNorm] = gcv(u,s,data,'Tikh', nLambda);
[lambdaO, gcvErrorMinNormO, lambdaGridMinNormO] = gcv(uO,sO,data,'Tikh', nLambda);

% Compute the minimum-norm estimates based on the gcv error.
betaMinNorm = cell(1, numel(lambdaGridMinNorm));
betaMinNormO = cell(1, numel(lambdaGridMinNorm));

for i = 1:numel(lambdaGridMinNorm)
    betaMinNorm{i} = zeros(size(templateOverlap, 2), size(data, 2));
    betaMinNormO{i} = zeros(size(templateOverlap, 2), size(data, 2));
    for ndx = 1:size(data, 2)
        betaMinNorm{i}(: ,ndx) = tikhonov( u , s , v , data(:, ndx), lambdaGridMinNorm(i));
        betaMinNormO{i}(: ,ndx) = tikhonov( uO , sO , vO , data(:, ndx), lambdaGridMinNormO(i));
    end
end

% Compute the betas for the "best" lambda from the gcv + other lambdas in
% the range
betaMinNormBest = zeros(size(betaMinNorm{1}));
betaMinNormAll = zeros([size(betaMinNorm{1}) nLambda]);
betaMinNormBestO = zeros(size(betaMinNorm{1}));
betaMinNormAllO = zeros([size(betaMinNorm{1}) nLambda]);
for i = 1:size(data, 2)
    betaMinNormBest(:, i) = tikhonov(u, s, v, data(:, i), lambda);
    betaMinNormBestO(:, i) = tikhonov(uO, sO, vO, data(:, i), lambdaO);
    for ll=1:nLambda
        betaMinNormAll(:, i,ll) = tikhonov(u, s, v, data(:, i), lambdaGridMinNorm(ll));
        betaMinNormAllO(:, i,ll) = tikhonov(uO, sO, vO, data(:, i), lambdaGridMinNormO(ll));
    end
end

% plot the gcv error 
figure;hold on;
plot(lambdaGridMinNorm,gcvErrorMinNorm,'LineWidth',2);xlabel('lambda');ylabel('gcv error');
plot(lambdaGridMinNormO,gcvErrorMinNormO,'LineWidth',2);xlabel('lambda');ylabel('gcv error');
legend('non-overlap','overlap','Location','Best')
saveas(gcf,'figures/gcvError.png','png')

% plot betas
figure; 
subplot(6,2,1);imagesc(betaMinNormBest);title(['no-overlap pickedLambda:' num2str(round(lambda))])
xlabel('time');  ylabel('ROI');
subplot(6,2,2);imagesc(betaMinNormBestO);title(['overlap pickedLambda:' num2str(round(lambdaO))])
aa=3;
for ll=[1 5 10 15 20]
    subplot(6,2,aa);
    imagesc(betaMinNormAll(:,:,ll));title(num2str(round(lambdaGridMinNorm(ll))))
    xlabel('time');  ylabel('ROI');
    subplot(6,2,aa+1);
    imagesc(betaMinNormAllO(:,:,ll));title(num2str(round(lambdaGridMinNormO(ll))))
    aa=aa+2;
end
set(gcf,'Position',[100 100 500 1000])
saveas(gcf,'figures/betasOverTime.png','png')


% %%%%%%%%%%%%%
% % compute fit and residuals (ind sbj) no scaling
% fitData = zeros([3200 90 5]); fitResidual = zeros([3200 90 5]);
% fitData(:,:,1) = templateNoOverlap*betaMinNormBest;
% fitResidual(:,:,1) = data-fitData(:,:,1);
% fitData(:,:,2) = templateNoOverlap*betaMinNormBestO;
% fitResidual(:,:,2) = data-fitData(:,:,2);
% pickL = [5 10 15];
% for ll = 1:3
%     fitData(:,:,ll+2) = templateNoOverlap*betaMinNormAll(:,:,pickL(ll));
%     fitResidual(:,:,ll+2) = data-fitData(:,:,ll+2);
% end
% % plot
% titleName = {'picked lambda','lambda from overlap set', ['lambda=' num2str(round(lambdaGridMinNorm(5)))],...
%     ['lambda=' num2str(round(lambdaGridMinNorm(10)))],['lambda=' num2str(round(lambdaGridMinNorm(15)))]};
% figure;
% subplot(6,2,1); plot(data(:,end));title('data');
% ff=3;
% for aa=1:5
%     subplot(6,2,ff); 
%     plot(squeeze(fitData(:,end,aa)));title(['fit ' titleName{aa}])
%     subplot(6,2,ff+1);
%     plot(squeeze(fitResidual(:,end,aa)));title('residuals')
%     ff=ff+2;
% end
% set(gcf,'Position',[100 100 1300 800])
% saveas(gcf,'figures/indFit.png','png')



%%%%%%%%%%%%%
% % compute fit and residuals (average) no scaling 
% test=reshape(data,[128,25,90]);
% dataAvg = squeeze(mean(test,2));
% % figure;plotOnEgi(dataAvg(:,end))
% fitAvg = zeros([128 90 5]); fitAvgRes = zeros([128 90 5]);
% fitAvg(:,:,1) = templateNoOverlap(1:128,:)*betaMinNormBest;
% fitAvgRes(:,:,1) = dataAvg-fitAvg(:,:,1);
% fitAvg(:,:,2) = templateNoOverlap(1:128,:)*betaMinNormBestO;
% fitAvgRes(:,:,2) = dataAvg-fitAvg(:,:,2);
% 
% fitAvgO = zeros([128 90 5]); fitResidualO = zeros([128 90 5]);
% fitAvgO(:,:,1) = templateOverlap(1:128,:)*betaMinNormBest;
% fitAvgResO(:,:,1) = dataAvg-fitAvgO(:,:,1);
% fitAvgO(:,:,2) = templateOverlap(1:128,:)*betaMinNormBestO;
% fitAvgResO(:,:,2) = dataAvg-fitAvgO(:,:,2);
% 
% pickL = [5 10 15];
% for ll = 1:3
%     fitAvg(:,:,ll+2) = templateNoOverlap(1:128,:)*betaMinNormAll(:,:,pickL(ll));
%     fitAvgRes(:,:,ll+2) = dataAvg-fitAvg(:,:,ll+2);
%     fitAvgO(:,:,ll+2) = templateOverlap(1:128,:)*betaMinNormAll(:,:,pickL(ll));
%     fitAvgResO(:,:,ll+2) = dataAvg-fitAvgO(:,:,ll+2);
% end
% 
% % plot
% titleName = {['bestLambdaNonOverlap=' num2str(round(lambda))],['bestLambdaOverlap=' num2str(round(lambdaO))], ...
%     ['lambda=' num2str(round(lambdaGridMinNorm(5)))],...
%     ['lambda=' num2str(round(lambdaGridMinNorm(10)))],['lambda=' num2str(round(lambdaGridMinNorm(15)))]};
% figure;
% subplot(6,2,1); plot(dataAvg(:,end),'LineWidth',2);title('average data');
% ff=3;
% for aa=1:5
%     subplot(6,2,ff); hold on
%     plot(squeeze(fitAvg(:,end,aa)),'LineWidth',2);
%     plot(squeeze(fitAvgO(:,end,aa)),'LineWidth',2);title(['fit ' titleName{aa}])
%     subplot(6,2,ff+1); hold on
%     plot(squeeze(fitAvgRes(:,end,aa)),'LineWidth',2);
%     plot(squeeze(fitAvgResO(:,end,aa)),'LineWidth',2);title('residuals')
%     ff=ff+2;
% end
% legend('non-overlap','overlap')
% set(gcf,'Position',[100 100 900 800])
% saveas(gcf,'figures/avgFit.png','png')
% 
% % plot just the residuals
% figure;
% subplot(2,1,1);hold on
% for aa=2:5
%     plot(squeeze(fitAvgRes(:,end,aa)),'LineWidth',2);
% end
% title('residuals non-overlap set')
% xlabel('electrodes')
% subplot(2,1,2);hold on
% for aa=2:5
%     plot(squeeze(fitAvgResO(:,end,aa)),'LineWidth',2);
% end
% legend(['lambda=' num2str(round(lambdaO))],['lambda=' num2str(round(lambdaGridMinNorm(5)))],...
%     ['lambda=' num2str(round(lambdaGridMinNorm(10)))],['lambda=' num2str(round(lambdaGridMinNorm(15)))]);
% xlabel('electrodes')
% title('residuals overlap set')
% saveas(gcf,'figures/residuals.png','png')









%%%%%%%
%%%%%%%%%%%%%
% compute fit and residuals (ind sbj) with scaling
fitData = zeros([3200 90 5]); fitResidual = zeros([3200 90 5]);
tempFit = templateNoOverlap*betaMinNormBest;
scalingCoef(1) = tempFit(:)\data(:);
fitData(:,:,1) = templateNoOverlap*(betaMinNormBest*scalingCoef(1));
fitResidual(:,:,1) = data-fitData(:,:,1);

tempFit = templateNoOverlap*betaMinNormBestO;
scalingCoef(2) = tempFit(:)\data(:);
fitData(:,:,2) = templateNoOverlap*(betaMinNormBestO*scalingCoef(2));
fitResidual(:,:,2) = data-fitData(:,:,2);

pickL = [5 10 15];
for ll = 1:3
    tempFit = templateNoOverlap*betaMinNormAll(:,:,pickL(ll));
    scalingCoef(ll+2) = tempFit(:)\data(:);
    fitData(:,:,ll+2) = templateNoOverlap*(betaMinNormAll(:,:,pickL(ll))*scalingCoef(ll+2));
    fitResidual(:,:,ll+2) = data-fitData(:,:,ll+2);
end
% plot
titleName = {'picked lambda','lambda from overlap set', ['lambda=' num2str(round(lambdaGridMinNorm(5)))],...
    ['lambda=' num2str(round(lambdaGridMinNorm(10)))],['lambda=' num2str(round(lambdaGridMinNorm(15)))]};
figure;
subplot(6,2,1); plot(data(:,end));title('data');
ff=3;
for aa=1:5
    subplot(6,2,ff); 
    plot(squeeze(fitData(:,end,aa)));title(['fit ' titleName{aa}])
    subplot(6,2,ff+1);
    plot(squeeze(fitResidual(:,end,aa)));title('residuals')
    ff=ff+2;
end
set(gcf,'Position',[100 100 1300 800])
saveas(gcf,'figures/indFit.png','png')




%%%%%%%%%%%%%
% compute fit and residuals (average) no scaling 
indData=reshape(data,[128,25,90]);
dataAvg = squeeze(mean(indData,2));
% figure;plotOnEgi(dataAvg(:,end))
fitAvg = zeros([128 90 5]); fitAvgRes = zeros([128 90 5]);
fitAvg(:,:,1) = templateNoOverlap(1:128,:)*(betaMinNormBest*scalingCoef(1));
fitAvgRes(:,:,1) = dataAvg-fitAvg(:,:,1);
fitAvg(:,:,2) = templateNoOverlap(1:128,:)*(betaMinNormBestO*scalingCoef(2));
fitAvgRes(:,:,2) = dataAvg-fitAvg(:,:,2);

% scaling for overlap set
tempFit = templateOverlap*betaMinNormBest;
scalingCoefO(1) = tempFit(:)\data(:);
tempFit = templateOverlap*betaMinNormBestO;
scalingCoefO(2) = tempFit(:)\data(:);
for ll = 1:3
    tempFit = templateOverlap*betaMinNormAllO(:,:,pickL(ll));
    scalingCoefO(ll+2) = tempFit(:)\data(:);
end

fitAvgO = zeros([128 90 5]); fitResidualO = zeros([128 90 5]);
fitAvgO(:,:,1) = templateOverlap(1:128,:)*(betaMinNormBest*scalingCoefO(1));
fitAvgResO(:,:,1) = dataAvg-fitAvgO(:,:,1);
fitAvgO(:,:,2) = templateOverlap(1:128,:)*(betaMinNormBestO*scalingCoefO(2));
fitAvgResO(:,:,2) = dataAvg-fitAvgO(:,:,2);

pickL = [5 10 15];
for ll = 3:5
    fitAvg(:,:,ll) = templateNoOverlap(1:128,:)*(betaMinNormAll(:,:,pickL(ll-2))*scalingCoef(ll));
    fitAvgRes(:,:,ll) = dataAvg-fitAvg(:,:,ll);
    fitAvgO(:,:,ll) = templateOverlap(1:128,:)*(betaMinNormAll(:,:,pickL(ll-2))*scalingCoef(ll));
    fitAvgResO(:,:,ll) = dataAvg-fitAvgO(:,:,ll);
end

% plot
titleName = {['bestLambdaNonOverlap=' num2str(round(lambda))],['bestLambdaOverlap=' num2str(round(lambdaO))], ...
    ['lambda=' num2str(round(lambdaGridMinNorm(5)))],...
    ['lambda=' num2str(round(lambdaGridMinNorm(10)))],['lambda=' num2str(round(lambdaGridMinNorm(15)))]};
figure;
subplot(6,2,1); plot(dataAvg(:,end),'LineWidth',2);title('average data');
ff=3;
for aa=1:5
    subplot(6,2,ff); hold on
    plot(squeeze(fitAvg(:,end,aa)),'LineWidth',2);
    plot(squeeze(fitAvgO(:,end,aa)),'LineWidth',2);title(['fit ' titleName{aa}])
    subplot(6,2,ff+1); hold on
    plot(squeeze(fitAvgRes(:,end,aa)),'LineWidth',2);
    plot(squeeze(fitAvgResO(:,end,aa)),'LineWidth',2);title('residuals')
    ff=ff+2;
end
legend('non-overlap','overlap')
set(gcf,'Position',[100 100 900 800])
saveas(gcf,'figures/avgFit.png','png')

% plot just the residuals
figure;
subplot(2,1,1);hold on
for aa=2:5
    plot(squeeze(fitAvgRes(:,end,aa)),'LineWidth',2);
end
title('residuals non-overlap set')
xlabel('electrodes')
subplot(2,1,2);hold on
for aa=2:5
    plot(squeeze(fitAvgResO(:,end,aa)),'LineWidth',2);
end
legend(['lambda=' num2str(round(lambdaO))],['lambda=' num2str(round(lambdaGridMinNorm(5)))],...
    ['lambda=' num2str(round(lambdaGridMinNorm(10)))],['lambda=' num2str(round(lambdaGridMinNorm(15)))]);
xlabel('electrodes')
title('residuals overlap set')
saveas(gcf,'figures/residuals.png','png')



test = templateNoOverlap(1:128,:)\squeeze(indData(:,1,:));




%%%%%%%%%% try l_curve
[reg_corner,rho,eta,reg_param] = l_curve(u,s,data(:,90),'Tikh');
[reg_cornerO,rhoO,etaO,reg_paramO] = l_curve(uO,sO,data(:,90),'Tikh');

betaLcurve = tikhonov(u, s, v, data(:, 90), reg_corner);
betaLcurveO = tikhonov(uO, sO, vO, data(:, 90), reg_cornerO);
figure;hold on;
plot(betaLcurve,'LineWidth',2); plot(betaLcurveO,'LineWidth',2)
plot(betaMinNormBest(:,90),'LineWidth',2); plot(betaMinNormBestO(:,90),'LineWidth',2)
xlabel('ROI')
legend(['lcurveDiff reg=' num2str(round(reg_corner))],['lcurveSame reg=' num2str(round(reg_cornerO))],...
    ['gcvDiff reg=' num2str(round(lambda))],['gcvSame reg=' num2str(round(lambdaO))])
saveas(gcf,'figures/lcurve.png','png')


tempFitLcurve = templateNoOverlap*betaLcurve;
scalingCoefLcurve = tempFitLcurve(:)\data(:,90);
fitDataLcurve = templateNoOverlap*(betaLcurve*scalingCoefLcurve);
fitResidualLcurve = data(:,90)-fitDataLcurve;
fitDataLcurveNos = templateNoOverlap*betaLcurve;
fitResidualLcurveNos = data(:,90)-fitDataLcurveNos;

tempFitLcurveO = templateOverlap*betaLcurveO;
scalingCoefLcurveO = tempFitLcurveO(:)\data(:,90);
fitDataLcurveO = templateOverlap*(betaLcurveO*scalingCoefLcurveO);
fitResidualLcurveO = data(:,90)-fitDataLcurveO;
fitDataLcurveONos = templateOverlap*(betaLcurveO);
fitResidualLcurveONos = data(:,90)-fitDataLcurveONos;

figure;hold on;
plot(betaLcurve,'LineWidth',2); plot(betaLcurveO,'LineWidth',2)
plot(betaLcurve*scalingCoefLcurve,'LineWidth',2); plot(betaLcurveO*scalingCoefLcurveO,'LineWidth',2)
legend('non-overlap','overlap','non-overlap scaled','overlap scaled')
xlabel('ROI')
saveas(gcf,'figures/lcurveScaling.png','png')
figure;hold on;
plot(fitResidualLcurveNos);plot(fitResidualLcurve)
plot(fitResidualLcurveONos);plot(fitResidualLcurveO)
rms(fitResidualLcurveNos) % not scaled
rms(fitResidualLcurve) % scaled
rms(fitResidualLcurveONos)
rms(fitResidualLcurveO)

