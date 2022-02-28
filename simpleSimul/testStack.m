clearvars;close all;
% Compare stack vs ind then average
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

totBoot = 10;

aucSta = zeros(length(snrLevel),totBoot);
energySta = aucSta;
mseStaNorm = aucSta;
aucInd = aucSta;
energyInd = aucSta;
mseAveInd = aucSta;

for noise = 1:length(snrLevel)
    noiseLevel = snrLevel(noise);
    
for repBoot=1:totBoot
    fprintf('noise %d bootstrap %d\n',noiseLevel,repBoot)
    
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
stackedForwards = [];
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
    stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[idxROIfwd{iSub,:}]));
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
regionROI = zeros(numSubs,numROIs,length(srcERP)); regionSta = regionROI;
for iSub=1:numSubs
    % regular minimum_norm: on the 20484 indexes per sbj
    [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    % beta values are for the indexes, but I want it per ROI
    % get the number of indexes per ROI for this subj
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    % get the range
    range = [0 cumsum(rangeROI)]; % cumulative sum of elements
    % SUM (not average) the beta values per ROI (=across the indexes)
    regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
end
% average across subj
retrieveROI = squeeze(mean(regionROI,1));

%% stacked ROI
stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));
tempY = permute(Y_avg,[2 1 3]);
stackY = reshape(tempY,[size(Y_avg,2)*numSubs,size(Y_avg,3)]);
stackY = bsxfun(@minus,stackY, mean(stackY));
[betaSta, lambdaSta] = minNormFast(stackedForwards, stackY, nLambdaRidge);
prevRange = 0; % for counting from prev sbj
for iSub=1:numSubs
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    range = [0 cumsum(rangeROI)] + prevRange;
    regionSta(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaSta(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    prevRange = range(end);
end
retrieveSta = squeeze(mean(regionSta,1));


    
[aucInd(noise,repBoot), energyInd(noise,repBoot), mseIndNorm(noise,repBoot)] = ...
    computeMetrics(squeeze(retrieveROI(:,winERP)),srcERP(:,winERP));
[aucSta(noise,repBoot), energySta(noise,repBoot), mseStaNorm(noise,repBoot)] = ...
    computeMetrics(squeeze(retrieveSta(:,winERP)),srcERP(:,winERP));

end
end

figure;hold on
subplot(1,3,1);hold on;
errorbar(log(snrLevel),mean(aucInd,2),std(aucInd,0,2),'LineWidth',2)
errorbar(log(snrLevel),mean(aucSta,2),std(aucSta,0,2),'LineWidth',2)
xlabel('log(SNR)');ylim([0 1]);ylabel('AUC')
subplot(1,3,2);hold on;
errorbar(log(snrLevel),mean(energyInd,2),std(energyInd,0,2),'LineWidth',2)
errorbar(log(snrLevel),mean(energySta,2),std(energySta,0,2),'LineWidth',2)
ylim([0 1]);ylabel('Energy');
subplot(1,3,3);hold on;
errorbar(log(snrLevel),mean(mseIndNorm,2),std(mseIndNorm,0,2),'LineWidth',2)
errorbar(log(snrLevel),mean(mseStaNorm,2),std(mseStaNorm,0,2),'LineWidth',2)
ylabel('MSE');ylim([0 1])
legend('ind','sta')
saveas(gcf,['figures/testStackROI'],'png')

