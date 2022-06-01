%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute auc, mse, relative energy & plot them

clearvars;close all;
addpath(genpath('subfunctions/'))

% V1 + MT
% load simulation results
load('simulOutput/simulV1MToutput.mat')

winERP = simulERP(1,1,1).winERP;
SNRlevel = unique([simulERP.noise]);
nbModel = 8;
% initialise variables
aucAve = zeros(size(simulERP,1),size(simulERP,2),size(simulERP,3),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;

for model=1:nbModel
for repBoot=1:size(simulERP,1)
    for totSbj=1:size(simulERP,2)
        for level=1:size(simulERP,3)            
        [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
            computeMetrics(squeeze(simulERP(repBoot,totSbj,level).beta(model,:,winERP)),simulERP(repBoot,totSbj,level).srcERP(:,winERP));        
        end
    end
end
end

%%% plot metrics
modName = {'template','whole','ROI','Oracle','wholeLC','ROILC','OracleLC','templateBestRegul'};

nbModel = 8;
figure;hold on
for model=1:nbModel
for ss=1:size(simulERP,2)
    subplot(3,nbModel,model);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
end
end
% legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 2200 700])
% set(gcf,'position',[100 100 1100 700])
saveas(gcf,['figures' filesep 'metricsV1MT'],'png')
saveas(gcf,['figures' filesep 'metricsV1MT'],'fig')
print(gcf,['figures' filesep 'metricsV1MT'],'-depsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot BETAs for sbj=50, noise=10, bootstrap=1
% noise=3;sbj=size(simulERP,2);repBoot=1;
% listROIs = simulERP(1,1,1).listROIs;
% 
% figure; count=1;
% currBeta = simulERP(repBoot,sbj,noise).srcERP;
% for iRoi = 1:2:length(listROIs)
%     subplot(5,9,count);hold on
%     plot(currBeta(iRoi,:) / max(max(abs(currBeta))) ,'LineWidth',2);
% %     plot(currBeta(iRoi+1,:) / max(max(abs(currBeta))) ,'LineWidth',2);
%     tt = cell2mat(listROIs(iRoi));title(tt(1:end-2) ,'LineWidth',2)
%     ylim([-1 1]);count=count+1;
% end
% % legend('left','right','location','best')
% for model=1:4
%     currBeta = squeeze(simulERP(repBoot,sbj,noise).beta(model,:,:));
%     for iRoi = 1:2:length(listROIs)
%         % need to normalise the signal
%         subplot(5,9,count);hold on
%         plot(currBeta(iRoi,:) / max(max(abs(currBeta))) ,'LineWidth',2);
% %         plot(currBeta(iRoi+1,:) / max(max(abs(currBeta))) ,'LineWidth',2);
%         tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%         ylim([-1 1]);count=count+1;
% %         legend('left','right')
%     end
% end

currBeta = simulERP(repBoot,sbj,noise).srcERP;
beta1 = squeeze(simulERP(repBoot,sbj,noise).beta(1,:,:));
beta2 = squeeze(simulERP(repBoot,sbj,noise).beta(2,:,:));
beta3 = squeeze(simulERP(repBoot,sbj,noise).beta(3,:,:));
beta4 = squeeze(simulERP(repBoot,sbj,noise).beta(4,:,:));
beta5 = squeeze(simulERP(repBoot,sbj,noise).beta(8,:,:));
count=1;
figure
for iRoi = 1:2:length(listROIs)
    subplot(6,9,count);hold on
    plot(currBeta(iRoi,:) / max(max(abs(currBeta))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2) ,'LineWidth',2)
    subplot(6,9,count+9);hold on
    plot(beta1(iRoi,:) / max(max(abs(beta1))) ,'LineWidth',2);ylim([-1 1]);
    title('template')
    subplot(6,9,count+9*2);hold on
    plot(beta2(iRoi,:) / max(max(abs(beta2))) ,'LineWidth',2);ylim([-1 1]);
    title('whole')
    subplot(6,9,count+9*3);hold on
    plot(beta3(iRoi,:) / max(max(abs(beta3))) ,'LineWidth',2);ylim([-1 1]);
    title('ROI')
    subplot(6,9,count+9*4);hold on
    plot(beta4(iRoi,:) / max(max(abs(beta4))) ,'LineWidth',2);ylim([-1 1]);
    title('oracle')
    subplot(6,9,count+9*5);hold on
    plot(beta5(iRoi,:) / max(max(abs(beta5))) ,'LineWidth',2);ylim([-1 1]);
    title('best regul')
    count = count+1;
end
set(gcf,'position',[100,100,2000,900])
saveas(gcf,'figures/betaErp' ,'png')
saveas(gcf,'figures/betaErp' ,'fig')
print(gcf,'figures/betaErp' ,'-depsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;close all;
load('simulOutput/simulV2V4output.mat')
winERP = simulERP(1,1).winERP;
nbModel = 8;
SNRlevel = unique([simulERP.noise]);
% initialise variables
aucAve = zeros(size(simulERP,1),nbModel,size(simulERP,2));
energyAve = aucAve;
mseAveNorm = aucAve;
for model=1:nbModel
for repBoot=1:size(simulERP,1)
    for totSbj=1:size(simulERP,2)
        for level=1:size(simulERP,3)            
        [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
            computeMetrics(squeeze(simulERP(repBoot,totSbj,level).beta(model,:,winERP)),simulERP(repBoot,totSbj,level).srcERP(:,winERP));        
        end
    end
end
end
%%% plot metrics
modName = {'template','whole','ROI','Oracle','wholeLC','ROILC','OracleLC','templateBestRegul'};
nbModel = 8;
figure;hold on
for model=1:nbModel
for ss=1:size(simulERP,2)
    subplot(3,nbModel,model);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
end
end
% legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 2200 700])
% set(gcf,'position',[100 100 1200 800])
saveas(gcf,['figures' filesep 'metricsV2V4'],'png')
saveas(gcf,['figures' filesep 'metricsV2V4'],'fig')
print(gcf,['figures' filesep 'metricsV2V4'],'-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot BETAs for sbj=50, noise=200, bootstrap=1
noise=4;sbj=size(simulERP,2);repBoot=1;
listROIs = simulERP(1,1,1).listROIs;

for model=1:4
    currBeta = squeeze(simulERP(repBoot,sbj,noise).beta(model,:,:));
    count = 1;
    figure;set(gcf,'position',[100,100,800,1000])
    for iRoi = 1:2:length(listROIs)
        % need to normalise the signal
        subplot(3,3,count);hold on
        plot(currBeta(iRoi,:) / max(max(abs(currBeta))) ,'LineWidth',2);
        plot(currBeta(iRoi+1,:) / max(max(abs(currBeta))) ,'LineWidth',2);
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    legend('left','right')
    saveas(gcf,['figures/betaErpV2V4' modName{model} '_SNR200'],'png')
%     saveas(gcf,['figures/betaErpV2V4' modName{model}],'fig')
%     saveas(gcf,['figures/betaErpV2V4' modName{model}],'eps')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;close all;
load('simulOutput/simulStepBilat.mat')
% 20 sbj, SNR=10, 2:2:18 active ROIs (bilateral activation)
winERP = simulBilat(1,1).winERP;
nbModel = 6; % models: 'template','whole','ROI','Oracle',OracleLC,'templateBestRegul'
% initialise variables
aucAve = zeros(size(simulBilat,1),nbModel,size(simulBilat,2));
energyAve = aucAve;
mseAveNorm = aucAve;
for model=1:nbModel
    for repBoot=1:size(simulBilat,1)
    for totROI=1:size(simulBilat,2)
        [aucAve(repBoot,model,totROI), energyAve(repBoot,model,totROI),mseAveNorm(repBoot,model,totROI)] = ...
            computeMetrics(squeeze(simulBilat(repBoot,totROI).beta(model,:,winERP)),simulBilat(repBoot,totROI).srcERP(:,winERP));        
    end
    end
end

%%% test SNR energy
% initialise variables
snr = zeros(size(simulBilat,1),nbModel,size(simulBilat,2));
snrPrct = snr;
winBaseline = 1:(min(winERP)-1);
for model=1:nbModel
    for repBoot=1:size(simulBilat,1)
        for totROI=1:size(simulBilat,2)
            [snr(repBoot,model,totROI), snrPrct(repBoot,model,totROI)]= computeMetricsSNR_test...
                (squeeze(simulBilat(repBoot,totROI).beta(model,:,:)),...
                simulBilat(repBoot,totROI).srcERP,winBaseline);
        end
    end
end

figure;hold on;
for model=1:6
    plot(2:2:18,log(mean(squeeze(snr(:,model,:)))),'LineWidth',2)
end
xlabel('nb active ROIs');ylabel('log SNR (A/C)')
legend('average','whole','ROI','OracleGCV','OracleLcurve','averageOracle','location','best');
saveas(gcf,['figures' filesep 'testSNRbilat'],'png')

figure;hold on;
for model=1:6
    plot(2:2:18,mean(squeeze(snrPrct(:,model,:))),'LineWidth',2)
end
xlabel('nb active ROIs');ylabel('log source % (A/(A+C)')
legend('average','whole','ROI','OracleGCV','OracleLcurve','averageOracle','location','best');
saveas(gcf,['figures' filesep 'testSNR%bilat'],'png')

% figure;hold on
% for mm=1:nbModel
% errorbar(2:2:18,squeeze(mean(pctEnergy(:,mm,:))),squeeze(std(pctEnergy(:,mm,:),1)),'LineWidth',2)
% end

%%% plot metrics
% figure;
% for mm=1:nbModel
%     subplot(2,2,1);hold on;
%     errorbar(2:2:18,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2)
%     xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([-1.5 4.5])xlim([2 18]);
%     subplot(2,2,2);hold on;
%     errorbar(2:2:18,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2)
%     ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5])xlim([2 18]);
%     subplot(2,2,3);hold on;
%     errorbar(2:2:18,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2)
%     ylim([0 1]);xlim([-1.5 4.5])ylabel('Energy');xlim([2 18]);
%     subplot(2,2,4);hold on;
%     errorbar(2:2:18,squeeze(mean(pctEnergy(:,mm,:))),squeeze(std(pctEnergy(:,mm,:),1)),'LineWidth',2)
%     ylim([0 1]);xlim([-1.5 4.5])ylabel('Energy SNR');xlim([2 18]);
% end
% legend('average','whole','ROI','OracleGCV','OracleLcurve');
figure;
for mm=[1:4 6]
    subplot(1,3,1);hold on;
    errorbar(2:2:18,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2,'CapSize',0)
    xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([2 18]);axis square
    subplot(1,3,2);hold on;
    errorbar(2:2:18,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);ylabel('Energy');xlim([2 18]);axis square
    subplot(1,3,3);hold on;
    errorbar(2:2:18,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1]);xlim([2 18]);axis square
end
legend('average','whole','ROI','Oracle','template optimal');
set(gcf,'position',[100 100 1000 300])
saveas(gcf,['figures' filesep 'nbOfBilatROIsStep'],'png')
saveas(gcf,['figures' filesep 'nbOfBilatROIsStep'],'fig')
print(gcf,['figures' filesep 'nbOfBilatROIsStep'],'-depsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;close all;
% load('simulOutput/simulRandUni.mat')
load('simulOutput/simulStepUni.mat')
% 20 sbj, SNR=10, 2:2:18 active ROIs (bilateral activation)
winERP = simulUni(1,1).winERP;
nbModel = 4;
% initialise variables
aucAve = zeros(size(simulUni,1),nbModel,size(simulUni,2));
energyAve = aucAve;
mseAveNorm = aucAve;
energyAveU = aucAve;
mseAveNormU = aucAve;
aucAveU=aucAve;
for model=1:nbModel
    for repBoot=1:size(simulUni,1)
    for totROI=1:size(simulUni,2)
        [aucAve(repBoot,model,totROI), energyAve(repBoot,model,totROI),mseAveNorm(repBoot,model,totROI)] = ...
            computeMetrics(squeeze(simulUni(repBoot,totROI).beta(model,:,winERP)),simulUni(repBoot,totROI).srcERP(:,winERP));        
        [aucAveU(repBoot,model,totROI), energyAveU(repBoot,model,totROI),mseAveNormU(repBoot,model,totROI)] = ...
            computeMetrics(squeeze(simulUni(repBoot,totROI).betaUni(model,:,winERP)),simulUni(repBoot,totROI).srcERPUni(:,winERP));
        % could try having 0 in the other ROI
        tempUni = zeros(18,length(winERP));
        if mod(repBoot,2)
            tempUni([1:2:18],:) = squeeze(simulUni(repBoot,totROI).betaUni(model,:,winERP));
        else
            tempUni([2:2:18],:) = squeeze(simulUni(repBoot,totROI).betaUni(model,:,winERP));
        end
        [aucAveU2(repBoot,model,totROI), energyAveU2(repBoot,model,totROI),mseAveNormU2(repBoot,model,totROI)] = ...
            computeMetrics(tempUni,simulUni(repBoot,totROI).srcERP(:,winERP));  
    end
    end
end
%%% plot metrics
figure;
for mm=1:nbModel
    subplot(3,3,1);hold on;
    errorbar(1:9,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2)
    xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([1 9]);
    title('18 potential ROI sources')
    subplot(3,3,2);hold on;
    errorbar(1:9,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');xlim([1 9]);
    subplot(3,3,3);hold on;
    errorbar(1:9,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1]);xlim([1 9]);
end
for mm=1:nbModel
    subplot(3,3,4);hold on;
    errorbar(1:9,squeeze(mean(aucAveU(:,mm,:))),squeeze(std(aucAveU(:,mm,:),1)),'LineWidth',2)
    xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([1 9]);
    title('9 potential ROI sources')
    subplot(3,3,5);hold on;
    errorbar(1:9,squeeze(mean(energyAveU(:,mm,:))),squeeze(std(energyAveU(:,mm,:),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');xlim([1 9]);
    subplot(3,3,6);hold on;
    errorbar(1:9,squeeze(mean(mseAveNormU(:,mm,:))),squeeze(std(mseAveNormU(:,mm,:),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1]);xlim([1 9]);
end
for mm=1:nbModel
    subplot(3,3,7);hold on;
    errorbar(1:9,squeeze(mean(aucAveU2(:,mm,:))),squeeze(std(aucAveU2(:,mm,:),1)),'LineWidth',2)
    xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([1 9]);
    title('9 ROI, metrics with 18 sources')
    subplot(3,3,8);hold on;
    errorbar(1:9,squeeze(mean(energyAveU2(:,mm,:))),squeeze(std(energyAveU2(:,mm,:),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');xlim([1 9]);
    subplot(3,3,9);hold on;
    errorbar(1:9,squeeze(mean(mseAveNormU2(:,mm,:))),squeeze(std(mseAveNormU2(:,mm,:),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1]);xlim([1 9]);
end
legend('average','whole','ROI','Oracle');
saveas(gcf,['figures' filesep 'nbOfUniROIsStep'],'png')



%%% plot BETAs for totROI=1 (2bilat)
repBoot=2; totROI = 1;
listROIs = simulUni(1,1,1).listROIs;
listROIsUni = listROIs(1:2:end);
tt = {'average','whole','ROI','Oracle'};
figure;
for model=1:4
    currBeta = squeeze(simulUni(repBoot,totROI).beta(model,:,min(winERP)));
    subplot(4,2,model+1*(model-1));hold on
    imagesc(currBeta)
    xticks(1:18);xticklabels(listROIs)
    title(tt(model))
    subplot(4,2,model+1*model);hold on
    currBeta = squeeze(simulUni(repBoot,totROI).betaUni(model,:,min(winERP)));
    imagesc(currBeta)
    xticks(1:9);xticklabels(listROIsUni)
end
saveas(gcf,['figures' filesep 'betaUniROIsStep'],'png')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;close all;

% load simulation results
load('simulOutput/simulV1MToutputStep.mat')

winERP = simulERP(1,1,1).winERP;
SNRlevel = unique([simulERP.noise]);
nbModel = 7;
% initialise variables
aucAve = zeros(size(simulERP,1),size(simulERP,2),size(simulERP,3),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;

for model=1:nbModel
for repBoot=1:size(simulERP,1)
    for totSbj=1:size(simulERP,2)
        for level=1:size(simulERP,3)            
        [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
            computeMetrics(squeeze(simulERP(repBoot,totSbj,level).beta(model,:,winERP)),simulERP(repBoot,totSbj,level).srcERP(:,winERP));        
        end
    end
end
end

%%% plot metrics
modName = {'average','whole','ROI','Oracle','wholeLC','ROILC','OracleLC'};

figure;hold on
for model=1:nbModel
for ss=1:size(simulERP,2)
    subplot(3,nbModel,model);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2)
    xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1])
end
end
legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 1500 700])
saveas(gcf,['figures' filesep 'metricsV1MTstep'],'png')













%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare systems
% compute auc, mse, relative energy & plot them
clearvars;close all;

% load simulation results
load('simulOutput/simulSysV1MT.mat')

winERP = simulSys(1,1,1).winERP;
SNRlevel = unique([simulSys.noise]);
nbModel = 5;

% initialise variables
aucAve = zeros(size(simulSys,1),size(simulSys,3),size(simulSys,2),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;

for model=1:nbModel
for sys=1:size(simulSys,3)
    for level=1:size(simulSys,2)
        for repBoot=1:size(simulSys,1)
            [aucAve(repBoot,level,sys,model), energyAve(repBoot,level,sys,model),mseAveNorm(repBoot,level,sys,model)] = ...
                computeMetrics(squeeze(simulSys(repBoot,level,sys).beta(model,:,winERP)),simulSys(repBoot,level,sys).srcERP(:,winERP));
        end
    end
end
end

%%% plot metrics
sysName = {'32','64','128','256'};

figure;hold on
for sys=1:size(simulSys,3)
for model=1:nbModel
    subplot(3,size(simulSys,3),sys);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,:,sys,model))),squeeze(std(aucAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
    xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
    title(sysName(sys))
    subplot(3,size(simulSys,3),sys+size(simulSys,3));hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,:,sys,model))),squeeze(std(energyAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(3,size(simulSys,3),sys+size(simulSys,3)*2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,:,sys,model))),squeeze(std(mseAveNorm(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1])
end
end
legend('Template','Whole','ROI','Oracle','TemplateBestRegul')
set(gcf,'position',[100 100 1500 700])
saveas(gcf,['figures' filesep 'compSysV1MT'],'png')

%%%%%%%%%%%%%%%%
% compare systems
% compute auc, mse, relative energy & plot them
clearvars;close all;

% load simulation results
load('simulOutput/simulSysV2V4.mat')

winERP = simulSys(1,1,1).winERP;
SNRlevel = unique([simulSys.noise]);
nbModel = 5;

% initialise variables
aucAve = zeros(size(simulSys,1),size(simulSys,3),size(simulSys,2),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;

for model=1:nbModel
for sys=1:size(simulSys,3)
    for level=1:size(simulSys,2)
        for repBoot=1:size(simulSys,1)
            [aucAve(repBoot,level,sys,model), energyAve(repBoot,level,sys,model),mseAveNorm(repBoot,level,sys,model)] = ...
                computeMetrics(squeeze(simulSys(repBoot,level,sys).beta(model,:,winERP)),simulSys(repBoot,level,sys).srcERP(:,winERP));
        end
    end
end
end

%%% plot metrics
sysName = {'32','64','128','256'};

figure;hold on
for sys=1:size(simulSys,3)
for model=1:nbModel
    subplot(3,size(simulSys,3),sys);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,:,sys,model))),squeeze(std(aucAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
    xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
    title(sysName(sys))
    subplot(3,size(simulSys,3),sys+size(simulSys,3));hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,:,sys,model))),squeeze(std(energyAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(3,size(simulSys,3),sys+size(simulSys,3)*2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,:,sys,model))),squeeze(std(mseAveNorm(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1])
end
end
legend('Template','Whole','ROI','Oracle','TemplateBestRegul')
set(gcf,'position',[100 100 1500 700])
saveas(gcf,['figures' filesep 'compSysV2V4'],'png')


%%%%%%%%%%%%%%%%
% compare systems
% compute auc, mse, relative energy & plot them
clearvars;close all;
% load simulation results
load('simulOutput/simulSysRand.mat')
winERP = simulSys(1,1,1).winERP;
SNRlevel = unique([simulSys.noise]);
nbModel = 5;
% initialise variables
aucAve = zeros(size(simulSys,1),size(simulSys,3),size(simulSys,2),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;
for model=1:nbModel
for sys=1:size(simulSys,3)
    for level=1:size(simulSys,2)
        for repBoot=1:size(simulSys,1)
            [aucAve(repBoot,level,sys,model), energyAve(repBoot,level,sys,model),mseAveNorm(repBoot,level,sys,model)] = ...
                computeMetrics(squeeze(simulSys(repBoot,level,sys).beta(model,:,winERP)),simulSys(repBoot,level,sys).srcERP(:,winERP));
        end
    end
end
end
%%% plot metrics for template
nameModel = {'Template','Whole','ROI','Oracle','TemplateBestRegul'};
nbModel = 5;
figure;hold on
for model=1:nbModel
    for sys=1:size(simulSys,3)
        subplot(nbModel,3,1+3*(model-1));hold on;
        errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,:,sys,model))),squeeze(std(aucAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
        xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC');xticks([-1:1:4]);axis square
        subplot(nbModel,3,2+3*(model-1));hold on;
        errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,:,sys,model))),squeeze(std(energyAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
        ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');xticks([-1:1:4]);axis square
        subplot(nbModel,3,3+3*(model-1));hold on;
        errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,:,sys,model))),squeeze(std(mseAveNorm(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
        ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);xticks([-1:1:4]);axis square
    end
%     legend({'32','64','128','256'})
    title(nameModel{model})
end
set(gcf,'position',[100 100 800 1000])
saveas(gcf,'figures/compSysRand','png')
saveas(gcf,'figures/compSysRand','fig')
print(gcf,'figures/compSysRand','-depsc')


%%%%% test SNR
load('simulOutput/simulV1MToutput.mat')

winERP = simulERP(1,1,1).winERP;
SNRlevel = unique([simulERP.noise]);
nbModel = 8;
% initialise variables
aucAve = zeros(size(simulERP,1),size(simulERP,2),size(simulERP,3),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;
winBaseline = 1:min(winERP)-1;

for model=1:nbModel
for repBoot=1:size(simulERP,1)
    for totSbj=1:size(simulERP,2)
        for level=1:size(simulERP,3)            
%         [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
%             computeMetrics(squeeze(simulERP(repBoot,totSbj,level).beta(model,:,winERP)),simulERP(repBoot,totSbj,level).srcERP(:,winERP));        
        [snr(repBoot,totSbj,level,model), snrPrct(repBoot,totSbj,level,model)]= computeMetricsSNR_test...
                (squeeze(simulERP(repBoot,totSbj,level).beta(model,:,:)),...
                simulERP(repBoot,totSbj,level).srcERP,winBaseline);        
        end
    end
end
end

figure;hold on;
for model=[1:4 8]
    plot(log([0.1 1 10 200 10000]),log(mean(squeeze(snr(:,size(simulERP,2),:,model)))),'LineWidth',2)
end
xlabel('log noise level');ylabel('log SNR (A/C')
legend({'average','whole','ROI','Oracle','averageOracle'},'location','best');
saveas(gcf,['figures' filesep 'testSNR'],'png')
figure;hold on;
for model=[1:4 8]
    plot(log([0.1 1 10 200 10000]),mean(squeeze(snrPrct(:,size(simulERP,2),:,model))),'LineWidth',2)
end
xlabel('log noise level');ylabel('source % (A/(A+C)')
legend({'average','whole','ROI','Oracle','averageOracle'},'location','best');
saveas(gcf,['figures' filesep 'testSNR%'],'png')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Cross Talk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load ('simulOutput/simulCrossTalk200.mat');
figure;
subplot(2,3,1);
imagesc(squeeze(mean(crossTalkTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))  
ylabel('seedArea');xlabel('predictArea')
title('templateBased')
subplot(2,3,2);imagesc(squeeze(mean(crossTalkWhole(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('Whole gcv')
subplot(2,3,3);imagesc(squeeze(mean(crossTalkROI(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('ROI gcv')
subplot(2,3,4);imagesc(squeeze(mean(crossTalkROIin(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('Oracle')
subplot(2,3,5);imagesc(squeeze(mean(crossTalkTemplateBest(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('template best reg param')
subplot(2,3,6);colorbar
colorcet('grey','reverse',1); % Gouldian reducedgrey heat L12 L18
colorcet('L12'); % Gouldian reducedgrey heat L12 L18
set(gcf,'position',[100,100,1500,1000])
saveas(gcf,['figures' filesep 'crossTalkStepSNR200' 'L12'],'png')
saveas(gcf,['figures' filesep 'crossTalkStepSNR200' 'L12'],'fig')
saveas(gcf,['figures' filesep 'crossTalkStepSNR200' 'L12'],'eps')
