clearvars;close all;
% unfinished program

addpath([pwd filesep 'subfunctions' filesep]);
load('averageMapEGI128.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
nLambdaRidge = 10;
lowPassNF1 = 1;

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
    Axx = load(['realData/eegdata/' sbjList(iSub).name '/Exp_MATL_HCN_128_Avg/Axx_c001.mat']);
    if lowPassNF1
        num2keep = 5;
        nf1 = Axx.i1F1;
        axxIdx = (nf1:nf1:10*nf1)+1;
        dftIdx = (1:10)+1;
        dft = dftmtx(size(Axx.Wave,1));
        sinWaveforms = imag(dft);
        cosWaveforms = real(dft);
        wave = cosWaveforms(:,dftIdx)*Axx.Cos(axxIdx,:)-sinWaveforms(:,dftIdx)*Axx.Sin(axxIdx,:);
        Axx.Wave = wave;
    end
    Y(iSub,:,:) = Axx.Wave';
end

%%% Use average reference for centering Y
%%% that is: substract the average electrode activity at each time point
% this is done by bsxfun which applies element-wise substraction (the 90
% averages across electrodes) - Useless
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end




%% compute minimum norm
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


regionWhole = zeros(numSubs,numROIs,length(Y_avg));
regionROI = zeros(numSubs,numROIs,length(Y_avg));
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
saveas(gcf,['figures/realTemplate' ],'png')
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
saveas(gcf,['figures/realROI' ],'png')
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
saveas(gcf,['figures/realWhole' ],'png')



count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(betaAverageModNorm(iRoi,:) / max(max(abs(betaAverageModNorm))) ,'LineWidth',2);
    plot(betaAverageModNorm(iRoi+1,:) / max(max(abs(betaAverageModNorm))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
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
saveas(gcf,['figures/realOracle' ],'png')
