clearvars;close all;
% Retrieve sources using real data from Lim 2017
% Bootstrap: rand sbj with replacement 500 times

% min norm on stalked sbj (same regularisation) is very long to repeat
% for all different set of participants for the bootstrap
% could do the min norm on separate sbj (different regularization) and then
% sample betas from that pool but regularisation won't be the same

addpath(genpath([pwd filesep 'subfunctions']))
load('template_EGI128.mat') % load average map of ROIs (128 elec x 18 ROIs)
avMap = templates.weights;
listROIs = templates.listROIs;
numROIs = length(listROIs);
nLambdaRidge = 50;
lowPassNF1 = 1; % filter 1 or not 0
numFq2keep = 5; % nb of harmonics to keep in the signal
totBoot = 500;
load('simulOutput2/realDataLasso.mat', 'sampleN')

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

%%% Re-ref data to average
for iSub=1:numSubs
    Y(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end


% bootstrap loop for 95% CI
betaTemplates = zeros(totBoot,numROIs,size(Y,3));
retrieveROI = betaTemplates;
retrieveWhole = betaTemplates;
retrieveROISum = betaTemplates;
retrieveWholeSum = betaTemplates;

for bb=1:totBoot
    clear regionROISum regionROI regionWholeSum regionWhole
    if mod(bb,10)==0;fprintf('boot nb %d',bb);end
%     pickN = randi(9,1,9);
    pickN = sampleN(bb,:);    
    sampleY = Y(pickN,:,:);
    
    %%%%%% "PCA"
    tempY = permute(sampleY,[2 1 3]);
    stackY = reshape(tempY,[size(Y,2)*numSubs,size(Y,3)]);
    stackY = bsxfun(@minus,stackY, mean(stackY));
    [u1, s1, v1] = svd(stackY);
    numComponents = 5;
    Ylo = u1(:,1:numComponents)*s1(1:numComponents,1:numComponents)*v1(:, 1:numComponents)';
    Y2 = reshape(Ylo,[size(Y,2),numSubs,size(Y,3)]);
    Ypca = permute(Y2,[2 1 3]);
    
    % get beta values for each ROI over time using templates
    betaTemplates(bb,:,:) = minNormFast_lcurve(avMap, squeeze(mean(Ypca,1)));
%     
%     %%%%%%%%%% Lim 2017 procedure
%     stackedFwd=[];stackedFullFwd = [];
%     for iSub=pickN
%         stackedFwd = blkdiag(stackedFwd, [roiFwd{iSub,:}]);
%         stackedFullFwd = blkdiag(stackedFullFwd, [fullFwd{iSub}]);
%     end
%     stackedFwd = bsxfun(@minus,stackedFwd, mean(stackedFwd));
%     stackedFullFwd = bsxfun(@minus,stackedFullFwd, mean(stackedFullFwd));
%     betaROI = minNormFast(stackedFwd,Ylo,50);
%     betaWhole = minNormFast(stackedFullFwd, Ylo,50);
%     
%     prevRange = 0; % for counting from prev sbj
%     nbS=1; % sbj count
%     for iSub=pickN
%         % get the number of indexes per ROI for this subj
%         rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
%         range = [0 cumsum(rangeROI)] + prevRange;
%         regionROISum(nbS,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%         regionROI(nbS,:,:) = cell2mat(arrayfun(@(x) mean(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%         % whole brain use idxROIfwd & add 20484 for each sbj to account
%         % for stacked data
%         regionWholeSum(nbS,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x} + (nbS-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
%         regionWhole(nbS,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole(idxROIfwd{iSub,x} + (nbS-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
%         % increment
%         prevRange = range(end);nbS = nbS+1;
%     end
%     retrieveROI(bb,:,:) = mean(regionROI); % for comparison with sum
%     retrieveROISum(bb,:,:) = mean(regionROISum); % = code
%     retrieveWhole(bb,:,:) = mean(regionWhole);
%     retrieveWholeSum(bb,:,:) = mean(regionWholeSum);
    
end
save('realDataOutput/realDataBootstrapTemplates2.mat','betaTemplates','sampleN')

% save('realDataOutput/realDataBootstrap100.mat','betaTemplates','retrieveROI','retrieveROISum','retrieveWhole','retrieveWholeSum','sampleN')



