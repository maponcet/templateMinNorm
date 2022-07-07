clearvars;close all;
% Retrieve sources using real data from Lim 2017
% No bootstrap, plot topographies at the end

addpath(genpath([pwd filesep 'subfunctions']))
load('template_EGI128.mat') % load average map of ROIs (128 elec x 18 ROIs)
avMap = templates.weights;
listROIs = templates.listROIs;
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

%%% Re-ref data to average
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

% get beta values for each ROI over time using templates
betaTemplates = minNormFast_lcurve(avMap, squeeze(mean(Ypca,1)));

%%%%%%%%%% Lim 2017 procedure
stackedFwd=[]; stackedFullFwd =[];
for iSub=1:numSubs
    stackedFwd = blkdiag(stackedFwd, [roiFwd{iSub,:}]);
    stackedFullFwd = blkdiag(stackedFullFwd, [fullFwd{iSub}]);
end
stackedFwd = bsxfun(@minus,stackedFwd, mean(stackedFwd));
[betaPlos,lambdaPlos] = minNormFast(stackedFwd,Ylo,50);
[betaWhole,lambdaWhole] = minNormFast(stackedFullFwd, Ylo,50);

prevRange = 0; % for counting from prev sbj
for iSub=1:numSubs
    % get the number of indexes per ROI for this subj
    rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
    range = [0 cumsum(rangeROI)] + prevRange;
    region(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    prevRange = range(end);
    % whole brain use idxROIfwd & add 20484 for each sbj to account
    % for stacked data
    regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole(idxROIfwd{iSub,x} + (iSub-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
end
retrieveROI = squeeze(mean(region)); 
retrieveWhole = squeeze(mean(regionWhole));

% reconstruct elec amplitude * time
reTemplates = avMap * betaTemplates; 
% for individual source loc, need to reconstruct per sbj
YhatPlos=stackedFwd*betaPlos;
unstackedYhatPlos = reshape(YhatPlos,128,numSubs,780);
reROI = squeeze(mean(unstackedYhatPlos,2));
YhatWhole=stackedFullFwd*betaWhole;
unstackedYhatWhole = reshape(YhatWhole,128,numSubs,780);
reWhole = squeeze(mean(unstackedYhatWhole,2));

save('realDataOutput/realData.mat','betaTemplates','retrieveROI','retrieveWhole','reTemplates','reROI','reWhole','YhatPlos','YhatWhole','Y','Ypca','t')
