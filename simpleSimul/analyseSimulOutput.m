%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute auc, mse, relative energy & plot them

clearvars;close all;
addpath(genpath('subfunctions/'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1 + MT
nbBoot = 30;
winERP = 46:180;
SNRlevel = [0.1 1 10 200 10000];
nbSbjToInclude =[2 8 20 50];
nbModel = 8;
% initialise variables
aucAve = zeros(nbBoot,length(nbSbjToInclude),length(SNRlevel),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;
for repBoot = 1:nbBoot
    % load simulation results
    fprintf('load bootstrap %d\n',repBoot)
    load(['simulOutput/simulV1MTfirst/simulV1MToutput' num2str(repBoot) '.mat'])
    for model=1:nbModel
        for totSbj=1:length(nbSbjToInclude)
            for level=1:length(SNRlevel)
                [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
                    computeMetrics(squeeze(simulERP(totSbj,level).beta(model,:,winERP)),simulERP(totSbj,level).srcERP(:,winERP));
            end
        end
    end
end

%%% plot metrics
modName = {'template','whole','ROI','Oracle','wholeLC','ROILC','OracleLC','templateBestRegul'};
figure;hold on
for model=1:nbModel
for ss=1:length(nbSbjToInclude)
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
legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 2200 700])
saveas(gcf,['figures' filesep 'metricsV1MT'],'png')
saveas(gcf,['figures' filesep 'metricsV1MT'],'fig')
print(gcf,['figures' filesep 'metricsV1MT'],'-depsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot BETAs for sbj=50, noise=10
sbj=length(nbSbjToInclude); noise=find(SNRlevel==10);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
currBeta = simulERP(sbj,noise).srcERP;
beta1 = squeeze(simulERP(sbj,noise).beta(1,:,:));
beta2 = squeeze(simulERP(sbj,noise).beta(2,:,:));
beta3 = squeeze(simulERP(sbj,noise).beta(3,:,:));
beta4 = squeeze(simulERP(sbj,noise).beta(4,:,:));
beta5 = squeeze(simulERP(sbj,noise).beta(8,:,:));
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
% V2 + V4
clearvars;close all;
nbBoot = 30;
winERP = 46:180;
SNRlevel = [0.1 1 10 200 10000];
nbSbjToInclude =[2 8 20 50];
nbModel = 8;
% initialise variables
aucAve = zeros(nbBoot,length(nbSbjToInclude),length(SNRlevel),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;
for repBoot = 1:nbBoot
    % load simulation results
    fprintf('load bootstrap %d\n',repBoot)
    load(['simulOutput/simulV2V4/simulV2V4output' num2str(repBoot) '.mat'])
    for model=1:nbModel
        for totSbj=1:length(nbSbjToInclude)
            for level=1:length(SNRlevel)
                [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
                    computeMetrics(squeeze(simulERP(totSbj,level).beta(model,:,winERP)),simulERP(totSbj,level).srcERP(:,winERP));
            end
        end
    end
end
%%% plot metrics
modName = {'template','whole','ROI','Oracle','wholeLC','ROILC','OracleLC','templateBestRegul'};
nbModel = 8;
figure;hold on
for model=1:nbModel
for ss=1:length(nbSbjToInclude)
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
legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 2200 700])
saveas(gcf,['figures' filesep 'metricsV2V4'],'png')
saveas(gcf,['figures' filesep 'metricsV2V4'],'fig')
print(gcf,['figures' filesep 'metricsV2V4'],'-depsc')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot BETAs for sbj=50, noise=200
noise=find(SNRlevel==200);sbj=length(nbSbjToInclude);
load('averageMap50Sum.mat') % load average map to get ROI names
for model=1:4
    currBeta = squeeze(simulERP(sbj,noise).beta(model,:,:));
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
%     saveas(gcf,['figures/betaErpV2V4' modName{model} '_SNR200'],'png')
%     saveas(gcf,['figures/betaErpV2V4' modName{model}],'fig')
%     saveas(gcf,['figures/betaErpV2V4' modName{model}],'eps')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% brain noise
clearvars; close all
nbBoot = 30;
winERP = 46:180; 
SNRlevel = [0.1 1 10 200 10000];
nbSbjToInclude =[2 8 20 50];
nbModel = 6;
% initialise variables
aucAve = zeros(nbBoot,length(nbSbjToInclude),length(SNRlevel),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;
for repBoot = 1:nbBoot
    % load simulation results
    fprintf('load bootstrap %d\n',repBoot)
    load(['simulOutput/brainNoise/simulV1MTbrainNoise' num2str(repBoot) '.mat'])
    for model=1:nbModel
        for totSbj=1:length(nbSbjToInclude)
            for level=1:length(SNRlevel)
                [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
                    computeMetrics(squeeze(simulERP(totSbj,level).beta(model,:,winERP)),simulERP(totSbj,level).srcERP(:,winERP));
            end
        end
    end
end

%%% plot metrics
modName = {'templateElecNoise','template','whole','ROI','Oracle','templateBestRegul'};
figure;hold on
for model=1:nbModel
for ss=1:length(nbSbjToInclude)
    subplot(3,nbModel,model);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC');
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
end
end
legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 1700 700])
saveas(gcf,['figures' filesep 'brainNoiseV1MT'],'png')
saveas(gcf,['figures' filesep 'brainNoiseV1MT'],'fig')
print(gcf,['figures' filesep 'brainNoiseV1MT'],'-depsc')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nb of active bilateral ROIs
clearvars; close all
snrLevel = [1 10 200];
nbBoot = 30;
modName = {'template','whole','ROI','Oracle','OracleLC','templateBestRegul'};
figure;
for noise=1:length(snrLevel)
    % initialise variables
    aucAve = zeros(nbBoot,length(modName),9);
    energyAve = aucAve;
    mseAveNorm = aucAve;
    load(['simulOutput/nbROI/simulStepBilatSNR' num2str(snrLevel(noise)) '.mat'])
    for model=1:length(modName)
        for repBoot=1:size(simulBilat,1)
            for totROI=1:size(simulBilat,2)
                [aucAve(repBoot,model,totROI), energyAve(repBoot,model,totROI),mseAveNorm(repBoot,model,totROI)] = ...
                    computeMetrics(squeeze(simulBilat(repBoot,totROI).beta(model,:,:)),simulBilat(repBoot,totROI).srcERP);
            end
        end
    end
    for mm=[1:4 6]
        subplot(3,3,1+3*(noise-1));hold on;
        errorbar(2:2:18,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2,'CapSize',0)
        xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([0 16]);xticks(0:2:16);axis square
        subplot(3,3,2+3*(noise-1));hold on;
        errorbar(2:2:18,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2,'CapSize',0)
        ylim([0 1]);ylabel('Energy');xlim([0 16]);xticks(0:2:16);axis square
        subplot(3,3,3+3*(noise-1));hold on;
        errorbar(2:2:18,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2,'CapSize',0)
        ylabel('MSE');ylim([0 1]);xlim([0 16]);xticks(0:2:16);axis square
    end
    title(['SNR' num2str(snrLevel(noise))])
end
legend(modName([1:4 6]))
set(gcf,'position',[100 100 1000 1000])
saveas(gcf,['figures' filesep 'nbOfBilatROIsStep' ],'png')
saveas(gcf,['figures' filesep 'nbOfBilatROIsStep' ],'fig')
print(gcf,['figures' filesep 'nbOfBilatROIsStepSNR' ],'-depsc')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nb of active bilateral ROIs ERP
clearvars; close all
snrLevel = [1 10];
nbBoot = 30;
winERP = 46:90;
modName = {'template','whole','ROI','Oracle','OracleLC','templateBestRegul'};
figure;
for noise=1:length(snrLevel)
    % initialise variables
    aucAve = zeros(nbBoot,length(modName),9);
    energyAve = aucAve;
    mseAveNorm = aucAve;
    load(['simulOutput/nbROI/simulERPBilat' num2str(snrLevel(noise)) '.mat'])
    for model=1:length(modName)
        for repBoot=1:size(simulBilat,1)
            for totROI=1:size(simulBilat,2)
                [aucAve(repBoot,model,totROI), energyAve(repBoot,model,totROI),mseAveNorm(repBoot,model,totROI)] = ...
                    computeMetrics(squeeze(simulBilat(repBoot,totROI).beta(model,:,winERP)),simulBilat(repBoot,totROI).srcERP(:,winERP));
            end
        end
    end
    for mm=[1:4 6]
        subplot(2,3,1+3*(noise-1));hold on;
        errorbar(2:2:18,squeeze(mean(aucAve(:,mm,:))),squeeze(std(aucAve(:,mm,:),1)),'LineWidth',2,'CapSize',0)
        xlabel('nb of active ROI');ylabel('AUC');ylim([0 1]);xlim([0 16]);xticks(0:2:16);axis square
        subplot(2,3,2+3*(noise-1));hold on;
        errorbar(2:2:18,squeeze(mean(energyAve(:,mm,:))),squeeze(std(energyAve(:,mm,:),1)),'LineWidth',2,'CapSize',0)
        ylim([0 1]);ylabel('Energy');xlim([0 16]);xticks(0:2:16);axis square
        subplot(2,3,3+3*(noise-1));hold on;
        errorbar(2:2:18,squeeze(mean(mseAveNorm(:,mm,:))),squeeze(std(mseAveNorm(:,mm,:),1)),'LineWidth',2,'CapSize',0)
        ylabel('MSE');ylim([0 1]);xlim([0 16]);xticks(0:2:16);axis square
    end
    title(['SNR' num2str(snrLevel(noise))])
end
legend(modName([1:4 6]))
set(gcf,'position',[100 100 1000 700])
saveas(gcf,['figures' filesep 'nbOfBilatROIsERP' ],'png')
saveas(gcf,['figures' filesep 'nbOfBilatROIsERP' ],'fig')
print(gcf,['figures' filesep 'nbOfBilatROIsERP' ],'-depsc')



















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare EEG systems with different nb of electrodes
clearvars;close all;

nbBoot = 30;
SNRlevel = [0.1 1 10 200 10000];
nbModel = 5;nbSys = 4;
% initialise variables
aucAve = zeros(nbBoot,length(SNRlevel),nbSys,nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;
for nn=1:nbBoot
% load simulation results
fprintf('load bootstrap %d\n',repBoot)
load(['simulOutput/compareSyst/simulSysRand' num2str(nn) '.mat'])
winERP = simulSys(1,1,1).winERP;
for model=1:nbModel
    for sys=1:nbSys
        for level=1:length(SNRlevel)
            [aucAve(nn,level,sys,model), energyAve(nn,level,sys,model),mseAveNorm(nn,level,sys,model)] = ...
                computeMetrics(squeeze(simulSys(level,sys).beta(model,:,winERP)),simulSys(level,sys).srcERP(:,winERP));
        end
    end
end
end
%%% plot metrics for template
nameModel = {'Template','Whole','ROI','Oracle','TemplateBestRegul'};
for model=1:nbModel
    figure;hold on
    for sys=1:nbSys
        subplot(1,3,1);hold on;
        errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,:,sys,model))),squeeze(std(aucAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
        xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC');xticks([-1:1:4]);axis square
        subplot(1,3,2);hold on;
        errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,:,sys,model))),squeeze(std(energyAve(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
        ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');xticks([-1:1:4]);axis square
        subplot(1,3,3);hold on;
        errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,:,sys,model))),squeeze(std(mseAveNorm(:,:,sys,model),1)),'LineWidth',2,'CapSize',0)
        ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);xticks([-1:1:4]);axis square
    end
    legend({'32','64','128','256'})
    title(nameModel{model})
    set(gcf,'position',[100 100 800 300])
%     saveas(gcf,['figures/compSysRand-' nameModel{model}],'png')
%     saveas(gcf,['figures/compSysRand-' nameModel{model}],'fig')
%     print(gcf,['figures/compSysRand-' nameModel{model}],'-depsc')
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Cross Talk
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for nn=[1 10 200]
load (['simulOutput/crossTalk/crossTalkN' num2str(nn) '.mat']); %1 10 200
figure;
subplot(2,3,1);
imagesc(squeeze(mean(crossTalkNormTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))  
ylabel('seedArea');xlabel('predictArea')
title('templateBased')
subplot(2,3,2);imagesc(squeeze(mean(crossTalkNormWhole(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('Whole gcv')
subplot(2,3,3);imagesc(squeeze(mean(crossTalkNormROI(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('ROI gcv')
subplot(2,3,4);imagesc(squeeze(mean(crossTalkNormROIin(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('Oracle')
subplot(2,3,5);imagesc(squeeze(mean(crossTalkNormTemplateBest(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('template best reg param')
subplot(2,3,6);colorbar
colorcet('grey','reverse',1); % Gouldian reducedgrey heat L12 L18
% colorcet('L12'); % Gouldian reducedgrey heat L12 L18
set(gcf,'position',[100,100,1500,1000])
saveas(gcf,['figures' filesep 'crossTalkN' num2str(nn) ],'png')
saveas(gcf,['figures' filesep 'crossTalkN' num2str(nn)],'fig')
saveas(gcf,['figures' filesep 'crossTalkN' num2str(nn)],'eps')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare step ERP
clearvars;close all
SNRlevel = [1 10 200];
figure;
for nn=1:3
    load (['simulOutput/crossTalk/crossTalkN' num2str(SNRlevel(nn)) '.mat']); %1 10 200
    subplot(2,3,nn);
    imagesc(squeeze(mean(crossTalkNormTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square;%caxis([0 1])
    set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))
    set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))
    ylabel('seedArea');xlabel('predictArea')
    title(['templateSNR' num2str(SNRlevel(nn))])
end
clear crossTalkNormTemplate
for nn=1:3
    load (['simulOutput/crossTalk/crossTalkERP' num2str(SNRlevel(nn)) '.mat']); %1 10 200
    subplot(2,3,nn+3);
    imagesc(squeeze(mean(crossTalkNormTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square;%caxis([0 1])
    set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))
    set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))
    ylabel('seedArea');xlabel('predictArea')
    title(['templateERP-SNR' num2str(SNRlevel(nn))])
end
colorcet('grey','reverse',1);
set(gcf,'position',[100,100,1500,1000])
saveas(gcf,['figures' filesep 'crossTalkComparison'],'png')










%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% PLOT output of simulation with non-visual interfering ROI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars; close all
nameOutput = {'simulV1MTinterference','simulV1MTinterferenceERP','simulVV2V4interferenceERP'};
nameOut = {'interfereV1MT','interfereV1MTerp','interfereV2V4erp'};
nbModel = 2;
totBoot = 60;
SNRlevel = [0.1 1 10 200 10000];
for nn=1:length(nameOutput)
    % initialise variables
    aucAve = zeros(totBoot,length(SNRlevel),nbModel);
    energyAve = aucAve;
    mseAveNorm = aucAve;
    clear extROI
    for repBoot=1:totBoot
        load(['simulOutput/interference/' nameOutput{nn} num2str(repBoot) '.mat'])
        winERP = simulERP(1).winERP;
        for model=1:nbModel
            for level=1:length(SNRlevel)
                [aucAve(repBoot,level,model), energyAve(repBoot,level,model),mseAveNorm(repBoot,level,model)] = ...
                    computeMetrics(squeeze(simulERP(level).beta(model,:,winERP)),simulERP(level).srcERP(:,winERP));
            end
        end
        extROI(repBoot) = simulERP(1).nonVisualROI;
    end
    %%% plot metrics
    modName = {'template','interference'};
    figure;hold on
    for model=1:nbModel
        subplot(1,3,1);hold on;
        errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,:,model))),squeeze(std(aucAve(:,:,model),1)),'LineWidth',2,'CapSize',0)
        xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
        subplot(1,3,2);hold on;
        errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,:,model))),squeeze(std(energyAve(:,:,model),1)),'LineWidth',2,'CapSize',0)
        ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
        subplot(1,3,3);hold on;
        errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,:,model))),squeeze(std(mseAveNorm(:,:,model),1)),'LineWidth',2,'CapSize',0)
        ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
    end
    legend('template','interference');
    set(gcf,'position',[100 100 700 300])
    saveas(gcf,['figures' filesep nameOut{nn}],'png')

%     % name of external source with resulting AUC
%     table(extROI', aucAve(:,1,2))
end







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% crosstalk with another active ROI
clearvars;close all
SNRlevel = [1 10 200];

for nn=1:length(SNRlevel)
    clear single singleV4
    figure;
    % single ROI
    load (['simulOutput/crossTalk/crossTalkN' num2str(SNRlevel(nn)) '.mat']); %1 10 200
    single = squeeze(mean(crossTalkNormTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18])));
    % 2 ROI
    load (['simulOutput/crossTalk/crossTalkV1N' num2str(SNRlevel(nn)) '.mat']); %1 10 200
    % plot
    subplot(2,2,1);
    imagesc(squeeze(mean(crossTalkNormTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square;%caxis([0 1])
    title('V1+')
    % predict sum for V1-L
    subplot(2,2,2);imagesc((single(1,:) + single) ./ max(single(1,:) + single));axis square
    title('predict V1+')
    % load 2 ROI with V4
    load (['simulOutput/crossTalk/crossTalkV4N' num2str(SNRlevel(nn)) '.mat']); %1 10 200
    singleV4 = squeeze(mean(crossTalkNormTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18])));
    subplot(2,2,3);imagesc(squeeze(mean(crossTalkNormTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));axis square;%caxis([0 1])
    title('V4+')
    % predict sum for V4-L
    subplot(2,2,4);    imagesc((singleV4(11,:) + singleV4) ./ max(singleV4(11,:) + singleV4));axis square
    title('predict V4+')
    colorcet('grey','reverse',1);
    saveas(gcf,['figures' filesep 'crossTalkAddSNR' num2str(SNRlevel(nn))],'png')
end
