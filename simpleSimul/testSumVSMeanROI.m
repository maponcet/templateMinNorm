clearvars;close all;
% Compare sum vs mean of ROI sources... Need to understand how it
% influences betas/sbj, average betas and the metrics
% use code to simulate V1-MT

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
noiseLevel = 10; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

numSubs = 20;

totBoot = 10;

for repBoot=1:totBoot

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

% list of random sbj with replacement
listSub = randi(length(dirList),numSubs,1);

%% LOAD FWD
fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
for iSub=1:numSubs
    clear fwdMatrix roiInfo
    % fwd file
    load([dataPath dirList(listSub(iSub)).name])
    fullFwd{iSub} = fwdMatrix;
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

%%% Simulate scalp activity (Y)
% use the generated sources to simulate scalp activity for each sbj
% (using individual fwd model)
Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));

for iSub=1:numSubs
    % initialise matrix of source activity
    sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
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
end

%%% Use average reference for centering Y
%%% that is: substract the average electrode activity at each time point
% this is done by bsxfun which applies element-wise substraction (the 90
% averages across electrodes) - useless?
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end

%% compute minimum norm
regionWhole = zeros(numSubs,numROIs,length(srcERP));
regionROI = regionWhole;
regionROILC = regionWhole;
regionWholeLC = regionWhole;
betaROIin = regionWhole;
betaROIinLC = regionWhole;

% min_norm on average data: get beta values for each ROI over time
[betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));

indFwdROI_noise=[roiFwd{iSub,:}];
indData_noise=squeeze(Y_avg(iSub,:,:));
indFwd_noise=fullFwd{iSub};

for iSub=1:numSubs
%     % regular minimum_norm: on the 20484 indexes per sbj
    [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%     [betaWholeLC, lambdaWholeLC] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
    
    [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    [betaROILC, lambdaROILC] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
    
    % beta values are for the indexes, but needs it per ROI
    % get the number of indexes per ROI for this subj
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    % get the range
    range = [0 cumsum(rangeROI)]; % cumulative sum of elements
    % SUM (not average) the beta values per ROI (=across the indexes)
    regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionROILC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionROI_mean(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionROILC_mean(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    
%     % need to find the indexes for whole brain -> use idxROIfwd
%     % (no need to get the range)
    regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%     regionWholeLC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
    regionWhole_mean(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%     regionWholeLC_mean(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%     
%     % feed ROI per sbj instead of mesh = "oracle" (best possible recovery)
%     sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
%     [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%     [betaROIinLC(iSub,:,:), lambdaGridMinNormROIinLC] = minNormFast_lcurve(sbjROI, squeeze(Y_avg(iSub,:,:)));
end
% average across subj
retrieveWhole = squeeze(mean(regionWhole,1));
retrieveROI = squeeze(mean(regionROI,1));
% retrieveROIin = squeeze(mean(betaROIin,1));
% retrieveWholeLC = squeeze(mean(regionWholeLC,1));
retrieveROILC = squeeze(mean(regionROILC,1));
% retrieveROIinLC = squeeze(mean(betaROIinLC,1));

retrieveWhole_mean = squeeze(mean(regionWhole_mean,1));
retrieveROI_mean = squeeze(mean(regionROI_mean,1));
% retrieveWholeLC_mean = squeeze(mean(regionWholeLC_mean,1));
retrieveROILC_mean = squeeze(mean(regionROILC_mean,1));


%%%% plot Betas
% count = 1;
% figure;set(gcf,'position',[100,100,800,1000])
% for iRoi = 1:2:length(listROIs)
%     % need to normalise the signal
%     subplot(3,3,count);hold on
%     plot(srcERP(iRoi,:) / max(max(abs(srcERP))) ,'LineWidth',2);
%     plot(retrieveWhole(iRoi,:) / max(max(abs(retrieveWhole))) ,'LineWidth',2);
%     plot(retrieveWhole_mean(iRoi+1,:) / max(max(abs(retrieveWhole_mean))) ,'LineWidth',2);
%     tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%     ylim([-1 1]);count=count+1;
% end
% legend('source','sum','mean');
if mod(repBoot,2) == 1
    count = 1;
    figure;set(gcf,'position',[100,100,800,1000])
    for iRoi = 1:2:length(listROIs)
        % need to normalise the signal
        subplot(3,3,count);hold on
        plot(srcERP(iRoi,:) / max(max(abs(srcERP))) ,'LineWidth',2);
        plot(retrieveROI(iRoi,:) / max(max(abs(retrieveROI))) ,'LineWidth',2);
        plot(retrieveROI_mean(iRoi,:) / max(max(abs(retrieveROI_mean))) ,'LineWidth',2);
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    legend('source','sum','mean');
end
if repBoot == 1
    count = 1;
    figure;set(gcf,'position',[100,100,800,1000])
    for iRoi = 1:2:length(listROIs)
        % need to normalise the signal
        subplot(3,3,count);hold on
        plot(srcERP(iRoi,:) / max(max(abs(srcERP))) ,'LineWidth',2);
        plot(retrieveROILC(iRoi,:) / max(max(abs(retrieveROILC))) ,'LineWidth',2);
        plot(retrieveROILC_mean(iRoi,:) / max(max(abs(retrieveROILC_mean))) ,'LineWidth',2);
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    legend('source','sum','mean');    
end
% count = 1;
% figure;set(gcf,'position',[100,100,800,1000])
% for iRoi = 1:2:length(listROIs)
%     % need to normalise the signal
%     subplot(3,3,count);hold on
%     plot(srcERP(iRoi,:)  ,'LineWidth',2);
%     plot(retrieveROI(iRoi,:)  ,'LineWidth',2);
%     plot(retrieveROI_mean(iRoi+1,:)  ,'LineWidth',2);
%     tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%     count=count+1;
% end
% legend('source','sum','mean');


[aucAve(repBoot), energyAve(repBoot), mseAveNorm(repBoot)] = ...
    computeMetrics(squeeze(retrieveROI(:,winERP)),srcERP(:,winERP));
[aucAve_mean(repBoot), energyAve_mean(repBoot), mseAveNorm_mean(repBoot)] = ...
    computeMetrics(squeeze(retrieveROI_mean(:,winERP)),srcERP(:,winERP));
[aucAveLC(repBoot), energyAveLC(repBoot), mseAveNormLC(repBoot)] = ...
    computeMetrics(squeeze(retrieveROILC(:,winERP)),srcERP(:,winERP));
[aucAveLC_mean(repBoot), energyAveLC_mean(repBoot), mseAveNormLC_mean(repBoot)] = ...
    computeMetrics(squeeze(retrieveROILC_mean(:,winERP)),srcERP(:,winERP));

[aucW(repBoot), energyW(repBoot), mseWNorm(repBoot)] = ...
    computeMetrics(squeeze(retrieveWhole(:,winERP)),srcERP(:,winERP));
[aucW_mean(repBoot), energyW_mean(repBoot), mseWNorm_mean(repBoot)] = ...
    computeMetrics(squeeze(retrieveWhole_mean(:,winERP)),srcERP(:,winERP));

end


% figure;hold on;
% bar(1:4, [mean(aucAve) mean(aucAve_mean) mean(aucAveLC) mean(aucAveLC_mean)]);
% figure;hold on;
% plot(aucAve,'LineWidth',2);
% plot(aucAve_mean,'LineWidth',2);
% plot(aucAveLC,'LineWidth',2);
% plot(aucAveLC_mean,'LineWidth',2);
% figure;hold on;
% plot(energyAve,'LineWidth',2);
% plot(energyAve_mean,'LineWidth',2);
% plot(energyAveLC,'LineWidth',2);
% plot(energyAveLC_mean,'LineWidth',2);
% figure;hold on;
% plot(mseAveNorm,'LineWidth',2);
% plot(mseAveNorm_mean,'LineWidth',2);
% plot(mseAveNormLC,'LineWidth',2);
% plot(mseAveNormLC_mean,'LineWidth',2);
% 
% figure;
% subplot(1,3,1);hold on;
% plot(aucAve,'LineWidth',2);
% plot(aucAve_mean,'LineWidth',2);
% ylabel('AUC');xlabel('simulation nb');
% subplot(1,3,2);hold on;
% plot(energyAve,'LineWidth',2);
% plot(energyAve_mean,'LineWidth',2);
% ylabel('Energy');xlabel('simulation nb');
% subplot(1,3,3);hold on;
% plot(mseAveNorm,'LineWidth',2);
% plot(mseAveNorm_mean,'LineWidth',2);
% ylabel('MSE');xlabel('simulation nb');
% legend('sum ROI','mean ROI')

figure;
subplot(1,3,1);hold on;
plot(aucW,'LineWidth',2);
plot(aucW_mean,'LineWidth',2);
ylabel('AUC');xlabel('simulation nb');
subplot(1,3,2);hold on;
plot(energyW,'LineWidth',2);
plot(energyW_mean,'LineWidth',2);
ylabel('Energy');xlabel('simulation nb');
subplot(1,3,3);hold on;
plot(mseWNorm,'LineWidth',2);
plot(mseWNorm_mean,'LineWidth',2);
ylabel('MSE');xlabel('simulation nb');
legend('sum Whole','mean Whole')


figure;
subplot(1,3,1);hold on;
plot(aucAve,'r','LineWidth',2);
plot(aucAve_mean,':r','LineWidth',2);
plot(aucW,'b','LineWidth',2);
plot(aucW_mean,':b','LineWidth',2);
ylabel('AUC');xlabel('simulation nb');
subplot(1,3,2);hold on;
plot(energyAve,'r','LineWidth',2);
plot(energyAve_mean,':r','LineWidth',2);
plot(energyW,'b','LineWidth',2);
plot(energyW_mean,':b','LineWidth',2);
ylabel('Energy');xlabel('simulation nb');
subplot(1,3,3);hold on;
plot(mseWNorm,'r','LineWidth',2);
plot(mseWNorm_mean,':r','LineWidth',2);
plot(mseAveNorm,'b','LineWidth',2);
plot(mseAveNorm_mean,':b','LineWidth',2);
ylabel('MSE');xlabel('simulation nb');
legend('sum ROI','mean ROI','sum Whole','mean Whole')

