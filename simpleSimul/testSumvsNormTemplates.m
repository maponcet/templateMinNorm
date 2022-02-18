clearvars;close all;
% Compare template using sum vs. normalised sum
% random 2 pairs of areas active (in both hemispheres, total=4sources)

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% test normalised map
regParam = sqrt(sum(avMap.^2,1));
avMapNorm = bsxfun(@rdivide,avMap,regParam);

% some parameters
snrLevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm

numSubs = 20;

totBoot = 30;

aucAve = zeros(length(snrLevel),totBoot);
energyAve = aucAve;
mseAveNorm = aucAve;
aucAveN = aucAve;
energyAveN = aucAve;
mseAveNormN = aucAve;

for noise = 1:length(snrLevel)
    noiseLevel = snrLevel(noise);
    
for repBoot=1:totBoot

    %% Simulate sources
    roiActive = randsample(1:2:numROIs,2);
    sourceL = listROIs(roiActive);
    sourceR = listROIs(roiActive+1);
    activeROIs = [sourceL,sourceR];
    % find the ROI index corresponding to the activeROIs
    ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));
    
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
% normalised template
[betaAverageNorm, lambdaNorm] = minNormFast_lcurve(avMapNorm, squeeze(mean(Y_avg,1)));

% count = 1;
% figure;set(gcf,'position',[100,100,800,1000])
% for iRoi = 1:2:length(listROIs)
%     % need to normalise the signal
%     subplot(3,3,count);hold on
%     plot(srcERP(iRoi,:) / max(max(abs(srcERP))) ,'LineWidth',2);
%     plot(betaAverageNorm(iRoi,:) / max(max(abs(betaAverageNorm))) ,'LineWidth',2);
%     plot(betaAverage(iRoi+1,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
%     tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%     ylim([-1 1]);count=count+1;
% end
% legend('source','avg','norm');

[aucAve(noise,repBoot), energyAve(noise,repBoot), mseAveNorm(noise,repBoot)] = ...
    computeMetrics(squeeze(betaAverage(:,winERP)),srcERP(:,winERP));
[aucAveN(noise,repBoot), energyAveN(noise,repBoot), mseAveNormN(noise,repBoot)] = ...
    computeMetrics(squeeze(betaAverageNorm(:,winERP)),srcERP(:,winERP));

end
end

figure;hold on
subplot(1,3,1);hold on;
errorbar(log(snrLevel),mean(aucAve,2),std(aucAve,0,2),'LineWidth',2)
errorbar(log(snrLevel),mean(aucAveN,2),std(aucAve,0,2),'LineWidth',2)
xlabel('log(SNR)');ylim([0 1]);ylabel('AUC')
subplot(1,3,2);hold on;
errorbar(log(snrLevel),mean(energyAve,2),std(energyAve,0,2),'LineWidth',2)
errorbar(log(snrLevel),mean(energyAveN,2),std(energyAveN,0,2),'LineWidth',2)
ylim([0 1]);ylabel('Energy');
subplot(1,3,3);hold on;
errorbar(log(snrLevel),mean(mseAveNorm,2),std(mseAveNorm,0,2),'LineWidth',2)
errorbar(log(snrLevel),mean(mseAveNormN,2),std(mseAveNormN,0,2),'LineWidth',2)
ylabel('MSE');ylim([0 1])
legend('sum','norm')

