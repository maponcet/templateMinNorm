clearvars;close all;
% unfinished program

addpath(genpath([pwd filesep 'subfunctions']))
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
nLambdaRidge = 20;
lowPassNF1 = 1; % filter or not
numFq2keep = 10; % nb of harmonics to keep in the signal

%% load forward & EEG data
sbjList = dir('realData/eegdata/skeri*');
numSubs = length(sbjList);

fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
for iSub=1:numSubs
    clear fwdMatrix roiInfo Axx
    % fwd file
    load(['realData/forwardSolutions/forwardAndRois-' sbjList(iSub).name '.mat'])
    fullFwd{iSub} = fwdMatrix;
    % go through each ROI and save the corresponding fwdMesh values
    % corresponding to the indexes of that ROI
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
        % save the index for each ROI
        idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
    end
    % eegFile
    % PlosOne uses condNmbr = 8;
    Axx = load(['realData/eegdata/' sbjList(iSub).name '/Exp_MATL_HCN_128_Avg/Axx_c008.mat']);
    if lowPassNF1
        nf1 = Axx.i1F1;
        axxIdx = (nf1:nf1:numFq2keep*nf1)+1;
        dftIdx = (1:numFq2keep)+1;
        dft = dftmtx(size(Axx.Wave,1));
        sinWaveforms = imag(dft);
        cosWaveforms = real(dft);
        wave = cosWaveforms(:,dftIdx)*Axx.Cos(axxIdx,:)-sinWaveforms(:,dftIdx)*Axx.Sin(axxIdx,:);
        Axx.Wave = wave;
    end
    Y(iSub,:,:) = Axx.Wave';
end
t = 0:Axx.dTms:(Axx.nT-1)*Axx.dTms;
figure;plot(t,squeeze(mean(Y))','k')
saveas(gcf,['figures/realData' num2str(numFq2keep)],'png')

%%% Use average reference for centering Y
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end




%% compute minimum norm
% min_norm on average data: get beta values for each ROI over time
[betaAverage, lambda] = minNormFast_lcurve(avMapNorm, squeeze(mean(Y_avg,1)));

regionWhole = zeros(numSubs,numROIs,length(Y_avg));
regionROI = zeros(numSubs,numROIs,length(Y_avg));
regionROILC = regionWhole;
regionWholeLC= regionWhole;
betaROIin= regionWhole;
betaROIinLC= regionWhole;

for iSub=1:numSubs
    % regular minimum_norm: on the 20484 indexes per sbj
    [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    
    [betaWholeLC,lambdaWholeLC] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
    [betaROILC, lambdaROILC] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));    
    % beta values are for the indexes, but I want it per ROI
    % get the number of indexes per ROI for this subj
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    % get the range
    range = [0 cumsum(rangeROI)]; % cumulative sum of elements
    % SUM (not average) the beta values per ROI (=across the indexes)
    regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionROILC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    
    % need to find the indexes for whole brain -> use idxROIfwd
    % (no need to get the range)
    regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
    regionWholeLC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
    
    % feed ROI per sbj instead of mesh
    sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
    [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    [betaROIinLC(iSub,:,:), lambdaGridMinNormROIinLC] = minNormFast_lcurve(sbjROI, squeeze(Y_avg(iSub,:,:)));
end
% average across subj
retrieveWhole = squeeze(mean(regionWhole,1));
retrieveROI = squeeze(mean(regionROI,1));
retrieveROIin = squeeze(mean(betaROIin,3));
retrieveROIinLC = squeeze(mean(betaROIinLC,3));
retrieveWholeLC = squeeze(mean(regionWholeLC,1));
retrieveROILC = squeeze(mean(regionROILC,1));        


%% Plots 
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
saveas(gcf,['figures/realTemplate' num2str(numFq2keep)],'png')
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
saveas(gcf,['figures/realROI' num2str(numFq2keep)],'png')
count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(retrieveROILC(iRoi,:) / max(max(abs(retrieveROILC))) ,'LineWidth',2);
    plot(retrieveROILC(iRoi+1,:) / max(max(abs(retrieveROILC))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/realROILC' num2str(numFq2keep)],'png')
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
saveas(gcf,['figures/realWhole' num2str(numFq2keep)],'png')
count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(retrieveWholeLC(iRoi,:) / max(max(abs(retrieveWholeLC))) ,'LineWidth',2);
    plot(retrieveWholeLC(iRoi+1,:) / max(max(abs(retrieveWholeLC))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/realWholeLC' num2str(numFq2keep)],'png')
count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(retrieveROIin(iRoi,:) / max(max(abs(retrieveROIin))) ,'LineWidth',2);
    plot(retrieveROIin(iRoi+1,:) / max(max(abs(retrieveROIin))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/realOracle' num2str(numFq2keep)],'png')

count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(retrieveROIinLC(iRoi,:) / max(max(abs(retrieveROIinLC))) ,'LineWidth',2);
    plot(retrieveROIinLC(iRoi+1,:) / max(max(abs(retrieveROIinLC))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/realOracleLC' num2str(numFq2keep)],'png')

