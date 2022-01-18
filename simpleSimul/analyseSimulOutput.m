%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute auc, mse, relative energy & plot them
% plot the betas
clearvars;close all;

% load simulation results
load('simulOutput/simulV1MToutput.mat')

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
    errorbar(log(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2)
    xlabel('log(SNR)');ylim([0 1]);ylabel('AUC')
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1])
end
end
legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 1500 700])
saveas(gcf,['figures' filesep 'metrics'],'png')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot BETAs for sbj=50, noise=10, bootstrap=1
noise=3;sbj=size(simulERP,2);repBoot=1;
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
    saveas(gcf,['figures/betaErp' modName{model}],'png')
    saveas(gcf,['figures/betaErp' modName{model}],'fig')
    saveas(gcf,['figures/betaErp' modName{model}],'eps')
end
currBeta = simulERP(repBoot,sbj,noise).srcERP;
count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:length(listROIs)
    subplot(3,3,count);hold on
    plot(currBeta(iRoi,:) / max(max(abs(currBeta))) ,'LineWidth',2);
    plot(currBeta(iRoi+1,:) / max(max(abs(currBeta))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2) ,'LineWidth',2)
    ylim([-1 1]);count=count+1;
end
legend('left','right')
saveas(gcf,'figures/betaErpSource','png')
saveas(gcf,'figures/betaErpSource','fig')

count = 1;
yData = simulERP(repBoot,totSbj,level).data;
yAvg = squeeze(mean(yData));
figure;plotOnEgi(yAvg(:,50))
figure;plotOnEgi(yAvg(:,120))
figure;plotOnEgi(yAvg(:,180))


currBeta = simulERP(repBoot,sbj,noise).srcERP;
beta1 = squeeze(simulERP(repBoot,sbj,noise).beta(1,:,:));
beta2 = squeeze(simulERP(repBoot,sbj,noise).beta(2,:,:));
beta3 = squeeze(simulERP(repBoot,sbj,noise).beta(3,:,:));
beta4 = squeeze(simulERP(repBoot,sbj,noise).beta(4,:,:));
count=1;
figure
for iRoi = 1:2:length(listROIs)
    subplot(5,9,count);hold on
    plot(currBeta(iRoi,:) / max(max(abs(currBeta))) ,'LineWidth',2);
    subplot(5,9,count+9);hold on
    plot(beta1(iRoi,:) / max(max(abs(beta1))) ,'LineWidth',2);ylim([-1 1]);
    subplot(5,9,count+9*2);hold on
    plot(beta2(iRoi,:) / max(max(abs(beta2))) ,'LineWidth',2);ylim([-1 1]);
    subplot(5,9,count+9*3);hold on
    plot(beta3(iRoi,:) / max(max(abs(beta3))) ,'LineWidth',2);ylim([-1 1]);
    subplot(5,9,count+9*4);hold on
    plot(beta4(iRoi,:) / max(max(abs(beta4))) ,'LineWidth',2);ylim([-1 1]);
    count = count+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;close all;
load('simulOutput/simulV2V4output.mat')
% 20 sbj, SNR=10, 2:2:18 active ROIs (bilateral activation)
winERP = simulERP(1,1).winERP;
nbModel = 7;
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
modName = {'average','whole','ROI','Oracle','wholeLC','ROILC','OracleLC'};
figure;hold on
for model=1:nbModel
for ss=1:size(simulERP,2)
    subplot(3,nbModel,model);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2)
    xlabel('log(SNR)');ylim([0 1]);ylabel('AUC')
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1])
end
end
legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 1500 700])
saveas(gcf,['figures' filesep 'metricsV2V4'],'png')

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
nbModel = 5;
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
%%% plot metrics
figure;
for mm=1:nbModel
    subplot(1,3,1);hold on;
    errorbar(2:2:18,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2)
    xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([2 18]);
    subplot(1,3,2);hold on;
    errorbar(2:2:18,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');xlim([2 18]);
    subplot(1,3,3);hold on;
    errorbar(2:2:18,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1]);xlim([2 18]);
end
legend('average','whole','ROI','OracleGCV','OracleLcurve');
saveas(gcf,['figures' filesep 'nbOfBilatROIsStep'],'png')



%%% test metrics SNR 

% initialise variables
aucAve = zeros(size(simulBilat,1),nbModel,size(simulBilat,2));
energyAve = aucAve;
mseAveNorm = aucAve;
for model=1:nbModel
    for repBoot=1:size(simulBilat,1)
    for totROI=1:size(simulBilat,2)
        [aucAve(repBoot,model,totROI), energyAve(repBoot,model,totROI),mseAveNorm(repBoot,model,totROI)] = ...
            computeMetricsSNR_test(squeeze(simulBilat(repBoot,totROI).beta(model,:,:)),simulBilat(repBoot,totROI).srcERP);        
    end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars;close all;
load('simulOutput/simulRandUni.mat')
% 20 sbj, SNR=10, 2:2:18 active ROIs (bilateral activation)
winERP = simulUni(1,1).winERP;
nbModel = 4;
% initialise variables
aucAve = zeros(size(simulUni,1),nbModel+nbModel,size(simulUni,2));
energyAve = aucAve;
mseAveNorm = aucAve;
for model=1:nbModel
    for repBoot=1:size(simulUni,1)
    for totROI=1:size(simulUni,2)
        [aucAve(repBoot,model,totROI), energyAve(repBoot,model,totROI),mseAveNorm(repBoot,model,totROI)] = ...
            computeMetrics(squeeze(simulUni(repBoot,totROI).beta(model,:,winERP)),simulUni(repBoot,totROI).srcERP(:,winERP));        
        [aucAve(repBoot,model+nbModel,totROI), energyAve(repBoot,model+nbModel,totROI),mseAveNorm(repBoot,model+nbModel,totROI)] = ...
            computeMetrics(squeeze(simulUni(repBoot,totROI).betaUni(model,:,winERP)),simulUni(repBoot,totROI).srcERPUni(:,winERP));  
    end
    end
end
%%% plot metrics
figure('position',[100 100 600 600]);
for mm=1:nbModel
    subplot(2,3,1);hold on;
    errorbar(1:9,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2)
    xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([1 9]);
    title('18 potential ROI sources')
    subplot(2,3,2);hold on;
    errorbar(1:9,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');xlim([1 9]);
    subplot(2,3,3);hold on;
    errorbar(1:9,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1]);xlim([1 9]);
end
for mm=nbModel+1:nbModel*2
    subplot(2,3,4);hold on;
    errorbar(1:9,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2)
    xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([1 9]);
    title('9 potential ROI sources')
    subplot(2,3,5);hold on;
    errorbar(1:9,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');xlim([1 9]);
    subplot(2,3,6);hold on;
    errorbar(1:9,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1]);xlim([1 9]);
end
legend('average','whole','ROI','Oracle');
saveas(gcf,['figures' filesep 'nbOfUniROIs'],'png')

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
    errorbar(log(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2)
    xlabel('log(SNR)');ylim([0 1]);ylabel('AUC')
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2)
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
load('simulOutput/simulSysV2V4.mat')

winERP = simulSys(1,1,1).winERP;
SNRlevel = unique([simulSys.noise]);
nbModel = 4;

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
    errorbar(log(SNRlevel),squeeze(mean(aucAve(:,:,sys,model))),squeeze(std(aucAve(:,:,sys,model),1)),'LineWidth',2)
    xlabel('log(SNR)');ylim([0 1]);ylabel('AUC')
    title(sysName(sys))
    subplot(3,size(simulSys,3),sys+size(simulSys,3));hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyAve(:,:,sys,model))),squeeze(std(energyAve(:,:,sys,model),1)),'LineWidth',2)
    ylim([0 1]);ylabel('Energy');
    subplot(3,size(simulSys,3),sys+size(simulSys,3)*2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,:,sys,model))),squeeze(std(mseAveNorm(:,:,sys,model),1)),'LineWidth',2)
    ylabel('MSE');ylim([0 1])
end
end
legend('Template','Whole','ROI','Oracle')
set(gcf,'position',[100 100 1500 700])
saveas(gcf,['figures' filesep 'compSysV2V4'],'png')

        