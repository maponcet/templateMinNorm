clearvars;close all;
% compare variability from sbj and from noise amount
% list of sbj changes across bootstrap but same across noise levels

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
% dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V2V-L','V4-L'};
sourceR = {'V2V-R','V4-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

nbSbjToInclude =[2 8 20 50];
% nbSbjToInclude =20;

totBoot = 10;

% initialise variables
aucAve = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyAve = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseAveNorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));

aucWhole = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
aucROI = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyWhole = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyROI = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseWholeNorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseROINorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));

aucAveModNorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyAveModNorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseAveNormModNorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
aucROIin = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyROIin = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseROINormin = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));

for totSbj=1:length(nbSbjToInclude)
numSubs = nbSbjToInclude(totSbj);


for repBoot=1:totBoot
    fprintf('N%d bootstrap %d\n',numSubs,repBoot)
% list of random sbj with replacement
listSub = randi(length(dirList),numSubs,1);
% since everything is fixed should sample without replacement -> RANDPERM
% (otherwise same sbj means that it's twice the same)

%% LOAD FWD
fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
for iSub=1:numSubs
    clear fwdMatrix roiInfo
    % fwd file
    load([dataPath dirList(listSub(iSub)).name])
    fullFwd{iSub} = fwdMatrix;
    %             indexROI = cell2mat(arrayfun(@(x) cellfind({roiInfo.name},listROIs{x}),1:length(listROIs),'uni',false));
    % go through each ROI and save the corresponding fwdMesh values
    % corresponding to the indexes of that ROI
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
        % to get roiFwd for one sbj= [roiFwd{iSub,:}]
        % save the index for each ROI
        idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
    end
end

%% Simulate sources
% amplitude (1 to 10) and time function is different for each
% source but the same for all sbj for a given bootstrap
% 1-45 = baseline 
% 46-90 = V1
% 91-135 = MT
% 136-180 = MT+V1
srcERP = zeros(numROIs,45*4); % 45*4 timepoints
srcERP(:,46:90) = createSourceERP(numROIs,ac_sources(1),ac_sources(3));
srcERP(:,91:135) = createSourceERP(numROIs,ac_sources(2),ac_sources(4));
srcERP(:,136:180) = createSourceERP(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));

% ERP & baseline timewindow
timeBase = 1:45;
winERP = 46:180;        


for level=1:length(SNRlevel)
% fprintf('noise %d\n',level)
%%
noiseLevel = SNRlevel(level);

%%% Simulate scalp activity (Y)
% use the generated sources to simulate scalp activity for each sbj
% (using individual fwd model)
Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
Y_noise = Y;
Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));

for iSub=1:numSubs
    % initialise matrix of source activity
    sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
    sourceNoise = sourceData; % no signal, used to compute SNR
    %             if length([idxROIfwd{iSub,ac_sources}]) ~= length(unique([idxROIfwd{iSub,ac_sources}]))
    %                 fprintf('Overlapping source indexes S%d \n',listSub(iSub));
    %             end
    for ss=1:length(ac_sources)
        % note that if there is overlapping index (same idx for 2
        % ROIs), the value in sourceData will be of the latest
        % source
        sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
    end
    % multiply fwd (128*20484) with the activated idx over time
    % (sourceData of 20484*90) and obtain Y elec x time
    y_stim = fullFwd{iSub} * sourceData;
    % add noise
    [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
    % to keep the same SNR for the 2 Y, need to compute noise for
    % the 2 Y separately as it is based on the variance of the signal
    Y(iSub,:,:) = y_stim + noisy_data;
    %             Y_SSVEPprev(iSub,:,:) = y_stimSSVEP + noisy_data;
    Y_noise(iSub,:,:) = noisy_data;
    %             Y_noiseSSVEP(iSub,:,:) = noisy_dataSSVEP;
    %             Ypure(iSub,:,:) = y_stim;
    %             Y2(iSub,:,:) = y_stim + noisy_dataSSVEP;
end

%%% Use average reference for centering Y
%%% that is: substract the average electrode activity at each time point
% this is done by bsxfun which applies element-wise substraction (the 90
% averages across electrodes) - Useless
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end


% figure;plotOnEgi(squeeze(mean(Y_avg(:,:,80),1)));colorbar
% figure;plotOnEgi(squeeze(mean(Y(:,:,80),1)));colorbar
% % % figure;plotOnEgi(squeeze(mean(Y_avgALL(:,:,80),1)));colorbar
% % figure;plotOnEgi(squeeze(Y_avg(1,:,80)));colorbar
% % figure;plotOnEgi(squeeze(Y(1,:,80)));colorbar
% timepoint = 50;
% figure;
% for iSub=1:numSubs
%     subplot(4,5,iSub)
%     plotOnEgi(Y_avg(iSub,:,timepoint)); colorbar
% end
% figure;
% for iSub=1:numSubs
%     subplot(4,5,iSub)
%     plotOnEgi(Y(iSub,:,timepoint)); colorbar
% end



%% compute minimum norm
regionWhole = zeros(numSubs,numROIs,length(srcERP));
regionROI = zeros(numSubs,numROIs,length(srcERP));
% min_norm on average data: get beta values for each ROI over time
% stack all the participants vertically
stackY = reshape(permute(Y_avg,[2,1,3]),[size(Y_avg,1)*size(Y_avg,2),size(Y_avg,3)]);
% stack the template for as many participants
stackAvMap = repmat(avMap,numSubs,1);
[betaAverage, betaMinNorm, lambda, gcvErrorMinNorm, lambdaGridMinNorm] = minimum_norm(stackAvMap, stackY, nLambdaRidge);


% do a normalisation for each ROI -> unit norming
% so that the betas represent microVolts (instead of microVolts/area size
% as it is now)
% unit norming is: all electrodes are squared and summed. These values are
% then divided so that the total of the electrodes for each ROI (power) is
% equal to 1
regParam = sqrt(sum(avMap.^2,1));
avMapNorm = bsxfun(@rdivide,avMap,regParam);
stackAvMapNorm = repmat(avMapNorm,numSubs,1);
[betaAverageModNorm, betaMinNormModNorm, lambdaModNorm, gcvErrorMinNormModNorm, lambdaGridMinNormModNorm] = ...
    minimum_norm(stackAvMapNorm, stackY, nLambdaRidge);



% for timepoint = 1:90
%     betaRegress(:,timepoint) = regress(squeeze(mean(Y_avg(:,:,timepoint),1))',avMap);
% end


% [betaAverage, ~, lambdaAverage(level)] = minimum_norm(avMap, squeeze(mean(Y_avg,1)), nLambdaRidge);
%         unstackedY = reshape(YloStack,size(fullFwd{1},1),numSubs,size(YloStack,2));
%         [betaAverageStack, ~, lambdaAverageStack] = minimum_norm(avMap, squeeze(mean(unstackedY,2)), nLambdaRidge);
%         unstackedSSVEP = reshape(Y_SSVEPloStack,size(fullFwd{1},1),numSubs,size(Y_SSVEPloStack,2));
%         [betaAverageSSVEP, ~, lambdaAverageSSVEP] = minimum_norm(avMap, squeeze(mean(unstackedSSVEP,2)), nLambdaRidge);
% figure;count = 1;
% for roi=1:2:18
%     subplot(3,3,count);hold on
%     plot(betaAverage(roi,:));
%     plot(betaAverageModNorm(roi,:));
%     count=count+1;
% end


% for iSub=1:numSubs
%     [betaAverage2(iSub,:,:), betaMinNorm2, lambda2(iSub), gcvErrorMinNorm2, lambdaGridMinNorm2] = minimum_norm(avMap, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
% end
% betaAverage2avg = squeeze(mean(betaAverage2,1));
% figure;count = 1;
% for roi=1:2:18
%     subplot(3,3,count);hold on
%     plot(betaAverage(roi,:));
%     plot(betaAverage2avg(roi,:));
%     count=count+1;
% end


for iSub=1:numSubs
    % regular minimum_norm: on the 20484 indexes per sbj
    [betaWhole, betaMinNormWhole, lambdaWhole,...
        gcvErrorMinNormWhole, lambdaGridMinNormWhole] = minimum_norm(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    % create brain resp 
    yWhole(iSub,:,:) = [fullFwd{iSub}] * betaWhole;

    [betaROI, betaMinNormROI, lambdaROI,...
        gcvErrorMinNormROI, lambdaGridMinNormROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    yROI(iSub,:,:) = [roiFwd{iSub,:}] * betaROI;

    % beta values are for the indexes, but I want it per ROI
    % get the number of indexes per ROI for this subj
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    % get the range
    range = [0 cumsum(rangeROI)]; % cumulative sum of elements
    % SUM (not average) the beta values per ROI (=across the indexes)
    regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');

    % need to find the indexes for whole brain -> use idxROIfwd
    % (no need to get the range)
    regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
    
    % test feeding ROI per sbj instead of mesh
    sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
    [betaROIin(:,:,iSub), betaMinNormROIin, lambdaROIin, gcvErrorMinNormROIin,...
        lambdaGridMinNormROIin] = minimum_norm(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
end
% average across subj
retrieveWhole = squeeze(mean(regionWhole,1));
retrieveROI = squeeze(mean(regionROI,1));
retrieveROIin = mean(betaROIin,3);



%% TOPO
%%%%%%%%%%%%%% 
% % plot topo for each sbj then plot the average of the topo
% timepoint = 50;
% figure;set(gcf,'position',[100,100,1500,1000])
% for iSub=1:numSubs
%     subplot(4,5,iSub)
%     plotOnEgi(yROI(iSub,:,timepoint)); colorbar
% end
% saveas(gcf,['figures/indTopoROIN' num2str(level)],'png')
% % figure;plotOnEgi(squeeze(mean(yROI(:,:,timepoint),1))); colorbar
% figure;set(gcf,'position',[100,100,1500,1000])
% for iSub=1:numSubs
%     subplot(4,5,iSub)
%     plotOnEgi(yWhole(iSub,:,timepoint)); colorbar
% end
% saveas(gcf,['figures/indTopoWholeN' num2str(level)],'png')
% % figure;plotOnEgi(squeeze(mean(yROI(:,:,timepoint),1))); colorbar
% % figure;plotOnEgi(squeeze(mean(yWhole(:,:,timepoint),1))); colorbar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute auc, mse, relative energy using average signal in rois
% do for all the min norm outputs
[aucAve(repBoot,totSbj,level), energyAve(repBoot,totSbj,level),...
    mseAveNorm(repBoot,totSbj,level),] = computeMetrics(betaAverage(:,winERP),srcERP(:,winERP));
[aucWhole(repBoot,totSbj,level), energyWhole(repBoot,totSbj,level),...
    mseWholeNorm(repBoot,totSbj,level)] = computeMetrics(retrieveWhole(:,winERP),srcERP(:,winERP));
[aucROI(repBoot,totSbj,level), energyROI(repBoot,totSbj,level),...
    mseROINorm(repBoot,totSbj,level)] = computeMetrics(retrieveROI(:,winERP),srcERP(:,winERP));

[aucAveModNorm(repBoot,totSbj,level), energyAveModNorm(repBoot,totSbj,level),...
    mseAveNormModNorm(repBoot,totSbj,level),] = computeMetrics(betaAverageModNorm(:,winERP),srcERP(:,winERP));
[aucROIin(repBoot,totSbj,level), energyROIin(repBoot,totSbj,level),...
    mseROINormin(repBoot,totSbj,level)] = computeMetrics(retrieveROIin(:,winERP),srcERP(:,winERP));

%% Plots for 1st bootstrap
if repBoot==1 && numSubs ==50
    %%% plot BETAs
    count = 1;
    figure;set(gcf,'position',[100,100,800,1000])
    for iRoi = 1:2:numROIs
        % need to normalise the signal
        subplot(3,3,count);hold on
        plot(srcERP(iRoi,:) / max(max(abs(srcERP))) ,'LineWidth',2);
        plot(srcERP(iRoi+1,:) / max(max(abs(srcERP))) ,'LineWidth',2);
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2) ,'LineWidth',2)
        ylim([-1 1]);count=count+1;
    end
    saveas(gcf,['figures/V2V4SrcN' num2str(numSubs) 'SNR' num2str(SNRlevel(level)) ],'png')
    count = 1;
    figure;set(gcf,'position',[100,100,800,1000])
    for iRoi = 1:2:numROIs
        % need to normalise the signal
        subplot(3,3,count);hold on
        plot(betaAverage(iRoi,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
        plot(betaAverage(iRoi+1,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    saveas(gcf,['figures/V2V4TempN' num2str(numSubs) 'SNR' num2str(SNRlevel(level)) ],'png')
    count = 1;
    figure;set(gcf,'position',[100,100,800,1000])
    for iRoi = 1:2:numROIs
        % need to normalise the signal
        subplot(3,3,count);hold on
        plot(retrieveROI(iRoi,:) / max(max(abs(retrieveROI))) ,'LineWidth',2);
        plot(retrieveROI(iRoi+1,:) / max(max(abs(retrieveROI))) ,'LineWidth',2);
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    saveas(gcf,['figures/V2V4RoiN' num2str(numSubs) 'SNR' num2str(SNRlevel(level)) ],'png')
    count = 1;
    figure;set(gcf,'position',[100,100,800,1000])
    for iRoi = 1:2:numROIs
        % need to normalise the signal
        subplot(3,3,count);hold on
        plot(retrieveWhole(iRoi,:) / max(max(abs(retrieveWhole))) ,'LineWidth',2);
        plot(retrieveWhole(iRoi+1,:) / max(max(abs(retrieveWhole))) ,'LineWidth',2);
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    saveas(gcf,['figures/V2V4WholeN' num2str(numSubs) 'SNR' num2str(SNRlevel(level)) ],'png')
    
%     %%% plot topo
%     % compare with average roiFwd for the participants then plot the average
%     % betas on this average roiFwd
%     roiMap = zeros(size(roiFwd{1},1),length(listROIs),numSubs);
%     for iSub=1:numSubs
%         for rr=1:length(listROIs)
%             roiMap(:,rr,iSub) = sum([roiFwd{iSub,rr}],2);
%         end
%     end
%     meanRoiMap = mean(roiMap,3);
%     
%     % plot topo
%     time1 = 70; time2=110; time3=160;
%     topoNew = avMap * betaAverage;
%     topoNewNorm = avMap * betaAverageModNorm;
%     topoWhole = meanRoiMap * retrieveWhole;
%     topoROI = meanRoiMap * retrieveROI;
%     figure;set(gcf,'position',[100,100,2300,1000])
%     subplot(3,8,1);plotOnEgi(squeeze(mean(Y(:,:,50),1)));% colorbar
%     subplot(3,8,2);plotOnEgi(squeeze(mean(Y_avg(:,:,50),1))); %colorbar
%     subplot(3,8,3);plotOnEgi(topoNew(:,50));% colorbar
%     subplot(3,8,4);plotOnEgi(topoWhole(:,50)); %colorbar
%     subplot(3,8,5);plotOnEgi(squeeze(mean(yWhole(:,:,50),1))); %colorbar
%     subplot(3,8,6);plotOnEgi(topoROI(:,50));% colorbar
%     subplot(3,8,7);plotOnEgi(squeeze(mean(yROI(:,:,50),1)));% colorbar
%     subplot(3,8,8);plotOnEgi(topoNewNorm(:,50));% colorbar
%     
%     subplot(3,8,9);plotOnEgi(squeeze(mean(Y(:,:,time2),1)));%colorbar
%     subplot(3,8,10);plotOnEgi(squeeze(mean(Y_avg(:,:,time2),1))); %colorbar
%     subplot(3,8,11);plotOnEgi(topoNew(:,time2)); %colorbar
%     subplot(3,8,12);plotOnEgi(topoWhole(:,time2)); %colorbar
%     subplot(3,8,13);plotOnEgi(squeeze(mean(yWhole(:,:,time2),1))); %colorbar
%     subplot(3,8,14);plotOnEgi(topoROI(:,time2)); %colorbar
%     subplot(3,8,15);plotOnEgi(squeeze(mean(yROI(:,:,time2),1))); %colorbar
%     subplot(3,8,16);plotOnEgi(topoNewNorm(:,time2)); %colorbar
%     
%     subplot(3,8,17);plotOnEgi(squeeze(mean(Y(:,:,time3),1))); %colorbar
%     subplot(3,8,18);plotOnEgi(squeeze(mean(Y_avg(:,:,time3),1))); %colorbar
%     subplot(3,8,19);plotOnEgi(topoNew(:,time3)); %colorbar
%     subplot(3,8,20);plotOnEgi(topoWhole(:,time3)); %colorbar
%     subplot(3,8,21);plotOnEgi(squeeze(mean(yWhole(:,:,time3),1))); %colorbar
%     subplot(3,8,22);plotOnEgi(topoROI(:,time3)); %colorbar
%     subplot(3,8,23);plotOnEgi(squeeze(mean(yROI(:,:,time3),1))); %colorbar
%     subplot(3,8,24);plotOnEgi(topoNewNorm(:,time3)); %colorbar
%     
%     saveas(gcf,['figures/topoErpN' num2str(numSubs) 'SNR' num2str(SNRlevel(level)) ],'png')
%     % order: 3 rows = V1, MT, V1+MT
%     % columns = source Y, average topo from each simul ind, template results,
%     % whole results from average betas, whole from average voltage (betas put
%     % into topo first for each sbj then average topo), roi results from betas,
%     % roi results from average voltages, normalised template
end

%%
end

end
end
% save('ERPtestSNR.mat','aucAve','energyAve','mseAveNorm','aucWhole','energyWhole','mseWholeNorm',...
%     'aucROI','energyROI','mseROINorm')


% %%% plot metrics
% figure;set(gcf,'position',[100,100,800,1200])
% subplot(3,3,1);hold on;
% plot(log(SNRlevel),squeeze(aucAve(:,1,:))','LineWidth',2)
% ylabel('AUC')
% ylim([0 1])
% title('template')
% subplot(3,3,2);hold on;
% plot(log(SNRlevel),squeeze(aucWhole(:,1,:))','LineWidth',2)
% ylim([0 1])
% title('whole')
% subplot(3,3,3);hold on;
% plot(log(SNRlevel),squeeze(aucROI(:,1,:))','LineWidth',2)
% ylim([0 1])
% title('ROI')
% 
% subplot(3,3,4);hold on;
% plot(log(SNRlevel),squeeze(energyAve(:,1,:))','LineWidth',2)
% ylabel('energy')
% ylim([0 1])
% subplot(3,3,5);hold on;
% plot(log(SNRlevel),squeeze(energyWhole(:,1,:))','LineWidth',2)
% ylim([0 1])
% subplot(3,3,6);hold on;
% plot(log(SNRlevel),squeeze(energyROI(:,1,:))','LineWidth',2)
% ylim([0 1])
% 
% subplot(3,3,7);hold on;
% plot(log(SNRlevel),squeeze(mseAveNorm(:,1,:))','LineWidth',2)
% ylabel('MSE')
% ylim([0 1])
% xlabel('log(SNR)')
% subplot(3,3,8);hold on;
% plot(log(SNRlevel),squeeze(mseWholeNorm(:,1,:))','LineWidth',2)
% ylim([0 1])
% subplot(3,3,9);hold on;
% plot(log(SNRlevel),squeeze(mseROINorm(:,1,:))','LineWidth',2)
% ylim([0 1])
% saveas(gcf,['figures' filesep 'bootVariabilityNtot'],'png')

%%% plot metrics
figure;
for ss=1:length(nbSbjToInclude)
    subplot(1,3,1);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucAve(:,ss,:))),squeeze(std(aucAve(:,ss,:),1)),'LineWidth',2)
    xlabel('log(SNR)')
    ylim([0 1])
    ylabel('AUC')
    title('template')
    subplot(1,3,2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyAve(:,ss,:))),squeeze(std(energyAve(:,ss,:),1)),'LineWidth',2)
    ylim([0 1])
    ylabel('Energy')
    subplot(1,3,3);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:))),squeeze(std(mseAveNorm(:,ss,:),1)),'LineWidth',2)
    ylabel('MSE')
    ylim([0 1])
end
legend('N2','N8','N20','N50')
saveas(gcf,['figures' filesep 'V2V4templateERP'],'png')

figure;
for ss=1:length(nbSbjToInclude)
    subplot(1,3,1);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucWhole(:,ss,:))),squeeze(std(aucWhole(:,ss,:),1)),'LineWidth',2)
    xlabel('log(SNR)')
    ylim([0 1])
    ylabel('AUC')
    title('Whole')
    subplot(1,3,2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyWhole(:,ss,:))),squeeze(std(energyWhole(:,ss,:),1)),'LineWidth',2)
    ylim([0 1])
    ylabel('Energy')
    subplot(1,3,3);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseWholeNorm(:,ss,:))),squeeze(std(mseWholeNorm(:,ss,:),1)),'LineWidth',2)
    ylabel('MSE')
    ylim([0 1])
end
legend('N2','N8','N20','N50')
saveas(gcf,['figures' filesep 'V2V4WholeERP'],'png')

figure;
for ss=1:length(nbSbjToInclude)
    subplot(1,3,1);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucROI(:,ss,:))),squeeze(std(aucROI(:,ss,:),1)),'LineWidth',2)
    xlabel('log(SNR)')
    ylim([0 1])
    ylabel('AUC')
    title('ROI')
    subplot(1,3,2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyROI(:,ss,:))),squeeze(std(energyROI(:,ss,:),1)),'LineWidth',2)
    ylim([0 1])
    ylabel('Energy')
    subplot(1,3,3);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseROINorm(:,ss,:))),squeeze(std(mseROINorm(:,ss,:),1)),'LineWidth',2)
    ylabel('MSE')
    ylim([0 1])
end
legend('N2','N8','N20','N50')
saveas(gcf,['figures' filesep 'V2V4NROIERP'],'png')


figure;
for ss=1:length(nbSbjToInclude)
    subplot(1,3,1);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucROIin(:,ss,:))),squeeze(std(aucROIin(:,ss,:))),'LineWidth',2)
    xlabel('log(SNR)')
    ylim([0 1])
    ylabel('AUC')
    title('ROIin')
    subplot(1,3,2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyROIin(:,ss,:))),squeeze(std(energyROIin(:,ss,:))),'LineWidth',2)
    ylim([0 1])
    ylabel('Energy')
    subplot(1,3,3);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseROINormin(:,ss,:))),squeeze(std(mseROINormin(:,ss,:))),'LineWidth',2)
    ylabel('MSE')
    ylim([0 1])
end
legend('N2','N8','N20','N50')
saveas(gcf,['figures' filesep 'V2V4ROIinERP'],'png')

% lineCOL={':r',':b',':g','--r','--b','--g','-r','-b','-g'};
% figure;
% subplot(1,3,1);hold on;
% for ss=1:length(nbSbjToInclude)
% errorbar(log(SNRlevel),squeeze(mean(aucAve(:,ss,:))),squeeze(std(aucAve(:,ss,:),1)),lineCOL{1+(ss-1)*3},'LineWidth',2)
% errorbar(log(SNRlevel),squeeze(mean(aucWhole(:,ss,:))),squeeze(std(aucWhole(:,ss,:))),lineCOL{2+(ss-1)*3},'LineWidth',2)
% errorbar(log(SNRlevel),squeeze(mean(aucROI(:,ss,:))),squeeze(std(aucROI(:,ss,:))),lineCOL{3+(ss-1)*3},'LineWidth',2)
% end
% xlabel('log(SNR)')
% ylabel('AUC')
% ylim([0 1])
% subplot(1,3,2);hold on;
% errorbar(log(SNRlevel),squeeze(mean(energyAve(:,ss,:))),squeeze(std(energyAve(:,ss,:))),lineCOL{1+(ss-1)*3},'LineWidth',2)
% errorbar(log(SNRlevel),squeeze(mean(energyWhole(:,ss,:))),squeeze(std(energyWhole(:,ss,:))),lineCOL{2+(ss-1)*3},'LineWidth',2)
% errorbar(log(SNRlevel),squeeze(mean(energyROI(:,ss,:))),squeeze(std(energyROI(:,ss,:))),lineCOL{3+(ss-1)*3},'LineWidth',2)
% xlabel('log(SNR)')
% ylabel('Energy')
% ylim([0 1])
% subplot(1,3,3);hold on;
% errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:))),squeeze(std(mseAveNorm(:,ss,:))),lineCOL{1+(ss-1)*3},'LineWidth',2)
% errorbar(log(SNRlevel),squeeze(mean(mseWholeNorm(:,ss,:))),squeeze(std(mseWholeNorm(:,ss,:))),lineCOL{2+(ss-1)*3},'LineWidth',2)
% errorbar(log(SNRlevel),squeeze(mean(mseROINorm(:,ss,:))),squeeze(std(mseROINorm(:,ss,:))),lineCOL{3+(ss-1)*3},'LineWidth',2)
% xlabel('log(SNR)')
% ylabel('MSEnorm')
% ylim([0 1])
% legend('20Template','20Whole','20ROI')
% set(gcf,'position',[100,100,1200,400])
% saveas(gcf,['figures' filesep 'SNRnTot'],'png')
    