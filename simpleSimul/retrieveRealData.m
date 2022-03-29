clearvars;close all;
% Retrieve sources using real data from Lim 2017

addpath(genpath([pwd filesep 'subfunctions']))
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
nLambdaRidge = 50;
lowPassNF1 = 1; % filter 1 or not 0
numFq2keep = 5; % nb of harmonics to keep in the signal

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


%% compute minimum norm

% min_norm on average data: get beta values for each ROI over time
[betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y,1)));
[betaAveragePCA, lambdaPCA] = minNormFast_lcurve(avMap, squeeze(mean(Ypca,1)));


%% Plots 
count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(t, betaAverage(iRoi,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
    plot(t, betaAverage(iRoi+1,:) / max(abs(betaAverage(:))) ,'LineWidth',2);
    line(t, zeros(size(t)),'Color','k')
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/realTemplate' num2str(numFq2keep)],'png')

count = 1;
figure;set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:numROIs
    % need to normalise the signal
    subplot(3,3,count);hold on
    plot(t, betaAveragePCA(iRoi,:) / max(max(abs(betaAveragePCA))) ,'LineWidth',2);
    plot(t, betaAveragePCA(iRoi+1,:) / max(abs(betaAveragePCA(:))) ,'LineWidth',2);
    line(t, zeros(size(t)),'Color','k')
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/realTemplatePCA' num2str(numFq2keep)],'png')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PlosOne procedure 
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
    plot(t, retrievePlosMean(iRoi,:) / max(max(abs(retrievePlosMean))) ,'LineWidth',2);
    plot(t, retrievePlosMean(iRoi+1,:) / max(max(abs(retrievePlosMean))) ,'LineWidth',2);
    line(t, zeros(size(t)),'Color','k')
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);count=count+1;
end
saveas(gcf,['figures/newPlos' num2str(numFq2keep)],'png')


