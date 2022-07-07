clearvars;close all;
% Test difference between 2 methods
% min norm on stalking sbj (same regularization) but is very long to repeat
% for all different set of participants for the bootstrap
% do the min norm on separate sbj (different regularization) and then
% sample betas from that pool = will be much much faster

addpath(genpath([pwd filesep 'subfunctions']))
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
% load('/Users/marleneponcet/Documents/Git/templateMinNorm/createTemplate/customEGI128.mat')
% avMap = customTemplates.data;
% listROIs = customTemplates.ROInames;
numROIs = length(listROIs);
nLambdaRidge = 50;
lowPassNF1 = 1; % filter 1 or not 0
numFq2keep = 5; % nb of harmonics to keep in the signal
totBoot = 500;

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
% figure;plot(t,squeeze(mean(Y))','k')
% saveas(gcf,['figures/realData' num2str(numFq2keep)],'png')

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


% compute min norm per sbj
regionROI = zeros(numROIs,length(Y),numSubs);
regionWhole = regionROI;
for iSub=1:numSubs
    % regular minimum_norm: on the 20484 indexes per sbj
    [betaWhole,lambdaWhole(iSub)] = minNormFast(fullFwd{iSub}, squeeze(Ypca(iSub,:,:)), nLambdaRidge);
    [betaROI, lambdaROI(iSub)] = minNormFast([roiFwd{iSub,:}], squeeze(Ypca(iSub,:,:)), nLambdaRidge);
    
    % beta values are for the indexes, but I want it per ROI
    % get the number of indexes per ROI for this subj
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    % get the range
    range = [0 cumsum(rangeROI)]; % cumulative sum of elements
    regionROI(:,:,iSub) = cell2mat(arrayfun(@(x) mean(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    
    % need to find the indexes for whole brain -> use idxROIfwd
    % (no need to get the range)
    regionWhole(:,:,iSub) = cell2mat(arrayfun(@(x) mean(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
end


% % bootstrap loop for 95% CI
% betaAverage = zeros(totBoot,numROIs,size(Y,3));
% betaAveragePCA = betaAverage;
% retrievePlos = betaAverage;
% retrievePlosMean = betaAverage;
% retrieveWhole = betaAverage;
betaROIresampling = zeros(totBoot,numROIs,size(Y,3));
betaWholeresampling = zeros(totBoot,numROIs,size(Y,3));
load('simulOutput2/realDataLasso.mat', 'sampleN')

for bb=1:totBoot
    clear region regionTmpMean regionWholeTmp
    if mod(bb,100)==0;fprintf('boot nb %d',bb);end
%     pickN = randi(9,1,9);
    pickN = sampleN(bb,:); 
    sampleY = Y(pickN,:,:);
    
    % quick way to resample for the min norm estimates but actually is
    % different if regularisation is done over the group (not on individual sbj) 
    betaROIresampling(bb,:,:) = mean(regionROI(:,:,pickN),3);
    betaWholeresampling(bb,:,:) = mean(regionWhole(:,:,pickN),3);
    
%     %%%%%% "PCA"
%     tempY = permute(sampleY,[2 1 3]);
%     stackY = reshape(tempY,[size(Y,2)*numSubs,size(Y,3)]);
%     stackY = bsxfun(@minus,stackY, mean(stackY));
%     [u1, s1, v1] = svd(stackY);
%     numComponents = 5;
%     Ylo = u1(:,1:numComponents)*s1(1:numComponents,1:numComponents)*v1(:, 1:numComponents)';
%     Y2 = reshape(Ylo,[size(Y,2),numSubs,size(Y,3)]);
%     Ypca = permute(Y2,[2 1 3]);
%     
%     
%     %% compute minimum norm
%     
%     % min_norm on average data: get beta values for each ROI over time
%     [betaAveragePCA(bb,:,:), lambdaPCA] = minNormFast_lcurve(avMap, squeeze(mean(Ypca,1)));
%     
%     %%%%%%%%%% PlosOne procedure
%     stackedFwd=[];stackedFullFwd = [];
%     for iSub=pickN
%         stackedFwd = blkdiag(stackedFwd, [roiFwd{iSub,:}]);
%         stackedFullFwd = blkdiag(stackedFullFwd, [fullFwd{iSub}]);
%     end
%     stackedFwd = bsxfun(@minus,stackedFwd, mean(stackedFwd));
%     stackedFullFwd = bsxfun(@minus,stackedFullFwd, mean(stackedFullFwd));
%     [betaPlos,lambdaPlosFull(bb)] = minNormFast(stackedFwd,Ylo,50);
%     [betaWhole,lambdaWholeFull(bb)] = minNormFast(stackedFullFwd, Ylo,50);
%     
%     prevRange = 0; % for counting from prev sbj
%     nbS=1; % sbj count
%     for iSub=pickN
%         % get the number of indexes per ROI for this subj
%         rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
%         range = [0 cumsum(rangeROI)] + prevRange;
%         region(nbS,:,:) = cell2mat(arrayfun(@(x) sum(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%         regionMean(nbS,:,:) = cell2mat(arrayfun(@(x) mean(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%         % whole brain use idxROIfwd & add 20484 for each sbj to account
%         % for stacked data
%         regionWholeTmp(nbS,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x} + (nbS-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
%         regionWholeTmpMean(nbS,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole(idxROIfwd{iSub,x} + (nbS-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
%         % increment
%         prevRange = range(end);nbS = nbS+1;
%     end
%     retrievePlos(bb,:,:) = mean(region); % for comparison with sum
%     retrievePlosMean(bb,:,:) = mean(regionMean); % = code
%     retrieveWhole(bb,:,:) = mean(regionWholeTmp);
%     retrieveWholeMean(bb,:,:) = mean(regionWholeTmpMean);
    
end

save('realDataOutput/realDataBootstrapResampling.mat','betaROIresampling','betaWholeresampling','sampleN')


%%%%% figures ROIs
% plot different outputs in separate line
count = 1;
plotMod = 2; figure;set(gcf,'position',[10,10,2400,600])
%     plotMod = 4; figure;set(gcf,'position',[10,10,2400,1000])
color = {'r','b'};
for oo=1:plotMod
    if oo==1
        data = betaWholeresampling;tname = 'whole';
    elseif oo==2
        data = betaROIresampling;tname = 'ROI';
    end
    % get the 95% CI
    ci95 = prctile(data,[2.5 97.5]);
    % value for normalising the signal
    normVal = max(abs(ci95(:)));
    % plot
    for iRoi = 1:numROIs
        subplot(plotMod,9,count);hold on
        plot(t, squeeze(mean(data(:,iRoi,:))) / normVal ,color{mod(iRoi,2)+1},'LineWidth',2);
        patch([t fliplr(t)], [squeeze(ci95(1,iRoi,:))' fliplr(squeeze(ci95(2,iRoi,:))')]/ normVal,...
            color{mod(iRoi,2)+1},'FaceAlpha',0.2, 'EdgeColor','none');
        line(t, zeros(size(t)),'Color','k')
        test0 = zeros(length(data),1);
        for tt=1:length(data)
            test0(tt) = sum(data(:,iRoi,tt)>0) / totBoot;
            if test0(tt) > 0.5
                test0(tt) = 1-test0(tt);
            end
        end
        sig0 = find((test0<0.025));
        if mod(iRoi,2)
            scatter(t(sig0), repmat(-0.8,length(sig0),1),'b.')
        else
            scatter(t(sig0), repmat(-0.9,length(sig0),1),'r.')
        end
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);
        if mod(iRoi,2) == 0
            count=count+1;
        end
    end
    title([tname tt(1:end-2)])
end
legend('left')
saveas(gcf,'figures/realDataResampling','png')
saveas(gcf,'figures/realDataResampling','fig')

