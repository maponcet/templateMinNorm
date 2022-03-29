%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Test stuff with real data
addpath(genpath([pwd filesep 'subfunctions']))
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
nLambdaRidge = 50;
lowPassNF1 = 1; % filter 1 or not 0
numFq2keep = 10; % nb of harmonics to keep in the signal

%%%%%%
%% LOAD AND FILTER THE EEG DATA (as specified) 
%%%%%%
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
        dftIdx = (1:numFq2keep)+1; % dftIdx has 1 step per Hz whereas axxIdx is 0.5Hz 
%         axxIdx = (nf1:nf1*2:numFq2keep*nf1)+1; 
%         dftIdx = (1:2:numFq2keep)+1; % dftIdx has 1 step per Hz whereas axxIdx is 0.5Hz 
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


%%% what reference is that data???
figure;hold on;
for iSub=1:numSubs
    plot(t,squeeze(mean(Y(iSub,:,:),2)))
end
% plot(t,squeeze(mean(mean(Y(:,:,:),2),1)),'k','linewidth',2)
title('mean electrodes for each sbj (after filtering non-harm)')
xlabel('time');
% saveas(gcf,['figures/realDataRefHarm'],'png')

%%% Re-ref to average
for iSub=1:numSubs
    Y(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end


%%%%%% "PCA"
tempY = permute(Y,[2 1 3]);
stackY = reshape(tempY,[size(Y,2)*numSubs,size(Y,3)]);
stackY = bsxfun(@minus,stackY, mean(stackY));
[u1, s1, v1] = svd(stackY);
numComponents = 5;
Ylo = u1(:,1:numComponents)*s1(1:numComponents,1:numComponents)*v1(:, 1:numComponents)';
Y2 = reshape(Ylo,[size(Y,2),numSubs,size(Y,3)]);
Ypca = permute(Y2,[2 1 3]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PlosOne procedure (be aware that minNormFast as a scalingCoef so
%%%%%%%%%% won't get the exact same betas)
stackedForwards=[];
for iSub=1:numSubs
    stackedForwards = blkdiag(stackedForwards, [roiFwd{iSub,:}]);
end
stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));
[betaPlos,lambdaPlos] = minNormFast(stackedForwards,Ylo,50);

prevRange = 0; % for counting from prev sbj
for iSub=1:numSubs
    % get the number of indexes per ROI for this subj 
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    range = [0 cumsum(rangeROI)] + prevRange;
    region(:,:,iSub) = cell2mat(arrayfun(@(x) sum(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionMean(:,:,iSub) = cell2mat(arrayfun(@(x) mean(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    prevRange = range(end);
end
retrievePlos = mean(region,3); % for comparison with sum 
retrievePlosMean = mean(regionMean,3); % = code

count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(retrievePlosMean(iRoi,:) / max(max(abs(retrievePlosMean))) ,'LineWidth',2);
    plot(retrievePlosMean(iRoi+1,:) / max(max(abs(retrievePlosMean))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/newPlos' num2str(numFq2keep)],'png')


% min_norm on average data: get beta values for each ROI over time
[betaAverage2, lambda2] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
[betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Ypca,1)));

regionWhole = zeros(numROIs,length(Ypca),numSubs);
regionROI = zeros(numROIs,length(Ypca),numSubs);
% regionROILC = regionWhole;
% regionWholeLC= regionWhole;
betaROIin= regionWhole;
% betaROIinLC= regionWhole;

% compared with stack data, lambda is different for each participant
for iSub=1:numSubs
    % regular minimum_norm: on the 20484 indexes per sbj
    [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Ypca(iSub,:,:)), nLambdaRidge);
    [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Ypca(iSub,:,:)), nLambdaRidge);
    
    [betaWholeLC,lambdaWholeLC] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
    [betaROILC, lambdaROILC] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));    
    % beta values are for the indexes, but I want it per ROI
    % get the number of indexes per ROI for this subj
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    % get the range
    range = [0 cumsum(rangeROI)]; % cumulative sum of elements
    regionROI(:,:,iSub) = cell2mat(arrayfun(@(x) mean(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionROILC(:,:,iSub) = cell2mat(arrayfun(@(x) mean(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    
    % need to find the indexes for whole brain -> use idxROIfwd
    % (no need to get the range)
    regionWhole(:,:,iSub) = cell2mat(arrayfun(@(x) mean(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
    regionWholeLC(:,:,iSub) = cell2mat(arrayfun(@(x) mean(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
    
    % feed ROI per sbj instead of mesh
    sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
    [betaROIin(:,:,iSub), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%     [betaROIinLC(iSub,:,:), lambdaGridMinNormROIinLC] = minNormFast_lcurve(sbjROI, squeeze(Y_avg(iSub,:,:)));
end
% average across subj
retrieveWhole = mean(regionWhole,3);
retrieveROI = mean(regionROI,3);
retrieveROIin = mean(betaROIin,3);
retrieveROILC = squeeze(mean(regionROILC,3));
retrieveWholeLC = squeeze(mean(regionWholeLC,3));
% retrieveROIinLC = squeeze(mean(regionROIinLC,1));        


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
% saveas(gcf,['figures/realTemplate_ICA_odd'],'png')
% saveas(gcf,['figures/realTemplate' num2str(numFq2keep)],'png')
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
% saveas(gcf,['figures/realROI_ICA_odd'],'png')

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
% saveas(gcf,['figures/realROILC_odd'],'png')

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
% saveas(gcf,['figures/realWhole_ICA_odd'],'png')

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
% saveas(gcf,['figures/realWholeLC_odd'],'png')

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
% saveas(gcf,['figures/realOracle_ICA_odd'],'png')
% saveas(gcf,['figures/realOracle' num2str(numFq2keep)],'png')

% count = 1;
% figure;set(gcf,'position',[100,100,800,1000])
% for iRoi = 1:2:numROIs
%     % need to normalise the signal
%     subplot(3,3,count);hold on
%     plot(retrieveROIinLC(iRoi,:) / max(max(abs(retrieveROIinLC))) ,'LineWidth',2);
%     plot(retrieveROIinLC(iRoi+1,:) / max(max(abs(retrieveROIinLC))) ,'LineWidth',2);
%     tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%     ylim([-1 1]);count=count+1;
% end
% saveas(gcf,['figures/realOracleLC' num2str(numFq2keep)],'png')





% compare different template methods
[betaSumICA, lambdaSumICA] = minNormFast_lcurve(avMap, squeeze(mean(Yica,1)));
[betaSumAvg, lambdaSumAvg] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));

meanAv = load('averageMap50Mean.mat') ;
avMapMean = meanAv.avMap.activity;
[betaMeanICA, lambdaMeanICA] = minNormFast_lcurve(avMapMean, squeeze(mean(Yica,1)));

% do a normalisation for each ROI -> unit norming
% so that the betas represent microVolts (instead of microVolts/area size
% as it is now)
% unit norming is: all electrodes are squared and summed. These values are
% then divided so that the total of the electrodes for each ROI (power) is
% equal to 1
regParam = sqrt(sum(avMap.^2,1));
avMapNorm = bsxfun(@rdivide,avMap,regParam);
regParamMean = sqrt(sum(avMapMean.^2,1));
avMapMeanNorm = bsxfun(@rdivide,avMapMean,regParamMean);
[betaSumICANorm, lambdaSumICANorm] = minNormFast_lcurve(avMapNorm, squeeze(mean(Yica,1)));
[betaMeanICANorm, lambdaMeanICANorm] = minNormFast_lcurve(avMapMeanNorm, squeeze(mean(Yica,1)));

% the 2 normalised maps (from sum or mean) are very very very similar
mm = round(max(max(abs(avMap))),-1);
figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMap(:,roi)) % only for 256 & 128 electrodes
    caxis([-mm mm])
end
saveas(gcf,['figures/avMap'],'png')

nn = round(max(max(abs(avMapNorm))),1);
figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMapNorm(:,roi)) % only for 256 & 128 electrodes
    caxis([-nn nn])
end
saveas(gcf,['figures/avMapNorm'],'png')

mm = round(max(max(abs(avMapMean))),-1);
figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMapMean(:,roi)) % only for 256 & 128 electrodes
    caxis([-mm mm])
end
saveas(gcf,['figures/avMapMean'],'png')

nn = round(max(max(abs(avMapMeanNorm))),1);
figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMapMeanNorm(:,roi)) % only for 256 & 128 electrodes
    caxis([-nn nn])
end
saveas(gcf,['figures/avMapMeanNorm'],'png')


plotBeta = betaMeanICANorm;
count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(plotBeta(iRoi,:) / max(max(abs(plotBeta))) ,'LineWidth',2);
    plot(plotBeta(iRoi+1,:) / max(max(abs(plotBeta))) ,'LineWidth',2);
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/realTestMeanICAnorm'],'png')


regParam = sqrt(sum(betaSumICA.^2,2));
betaSumICA_2 = bsxfun(@rdivide,betaSumICA,regParam);

