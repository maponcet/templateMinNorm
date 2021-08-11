clearvars;close all;
% simulate a pair of areas and compute how well the models do

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50.mat') % load average map of ROIs (128 elec x 18 ROIs)
listROIs = avMap.roiNames; % the order of the list MUST be the same as in the avMap!
numROIs = length(listROIs);
numSubs = length(dirList);

% some parameters
SNRlevel = 0.1; % noise level 10%
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm
numCols = 2; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
totBoot = 20; % nb of bootstrap

% allCombi = nchoosek(1:numROIs,2); % all possible combinations of areas

% initialise variables
aucNew = zeros(totBoot,length(listROIs),length(listROIs));
energyNew = zeros(totBoot,length(listROIs),length(listROIs));
mseNew = zeros(totBoot,length(listROIs),length(listROIs));
aucOld = zeros(totBoot,length(listROIs),length(listROIs));
energyOld = zeros(totBoot,length(listROIs),length(listROIs));
mseOld = zeros(totBoot,length(listROIs),length(listROIs));

for seed1 = 1:length(listROIs)
    for seed2 = 1:length(listROIs)
        if seed1==seed2
            break
        end
        
        signalROIs = [listROIs(seed1) listROIs(seed2)]; % simulated signal
        % find the ROI index which corresponds to the one in signalROIs
        ac_sources = cell2mat(arrayfun(@(x) cellfind(avMap.roiNames,signalROIs{x}),1:length(signalROIs),'uni',false));
        
        for repBoot =1:totBoot
            clear Y source signal Ylo
            fprintf('seed %d & seed %d bootstrap %d \n',seed1,seed2,repBoot);
            
            % list of random sbj with replacement
            listSub = randi(length(dirList),numSubs,1);
            
            %% GET ROIs and simulate signals
            % do for each sbj separately but on 18 ROIs, not mesh
            % so first need to get the ROI per sbj
            roiMap=zeros(128,numROIs,numSubs);
            roiStacked = [];
            for iSubj=1:numSubs
                clear fwdMatrix roiInfo
                load([dataPath dirList(listSub(iSubj)).name])
                % go through each ROI and average all the mesh indexes corresponding to that ROI
                for rr=1:numROIs
                    clear indexROI
                    indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
                    [roiMap(:,rr,iSubj)] = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
                end
                roiStacked = blkdiag(roiStacked, roiMap(:,:,iSubj));
            end
            roiStacked = bsxfun(@minus,roiStacked, mean(roiStacked));
            
            %% Simulate sources
           [Y, ~, sourceValOverTime, ~] = simulSource(roiMap,SNRlevel, ac_sources,[]);

            %%  reduce dimension data (Ylo)
            % =PCA denoised version of Y (denoised by truncation of the SVD)
            n = numel(Y);
            [u1, s1, v1] = svd(Y);
            Ylo = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
            % compute average EEG
            unstackedData = reshape(Ylo,128,numSubs,size(Ylo,2));
            grandMeanData = squeeze(mean(unstackedData,2));
            
            %% compute minimum norm
            [betaReg, ~, lambdaReg] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
            [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(roiStacked, Ylo, nLambdaRidge);
            % betaMinNorm is numSubs*numROIs x time so need to extract betaMinNorm of a single sbj
            % for each ROI separately then average across sbj
            regionTmp = arrayfun(@(x) betaMinNorm(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
            regionMinNorm =reshape(cell2mat(regionTmp),numROIs,size(Ylo,2),numSubs);
            regionMinNormAvg = mean(regionMinNorm,3);
            
            %% compute auc, mse, relative energy using average signal in rois
            [aucNew(repBoot,seed1,seed2), energyNew(repBoot,seed1,seed2),...
                mseNew(repBoot,seed1,seed2)] = computeMetrics(betaReg,ac_sources,sourceValOverTime);
%             [aucOld(repBoot,seed1,seed2), energyOld(repBoot,seed1,seed2),...
%                 mseOld(repBoot,seed1,seed2)] = computeMetrics(betaMinNorm,ac_sources,sourceValOverTime);
            [aucOld(repBoot,seed1,seed2), energyOld(repBoot,seed1,seed2),...
                mseOld(repBoot,seed1,seed2)] = computeMetrics(regionMinNormAvg,ac_sources,sourceValOverTime);
            
        end % end boostrap
        
    end % 
end
    
figure;
subplot(1,2,1);
imagesc(squeeze(mean(aucNew))); caxis([0 1]);colorbar
title('AUC New for retrieving pair of sources')
set(gca,'XTick',1:numROIs,'XTickLabel', listROIs );
set(gca,'YTick',1:numROIs,'YTickLabel', listROIs );
subplot(1,2,2);
imagesc(squeeze(mean(aucOld))); caxis([0 1]);colorbar
set(gca,'XTick',1:numROIs,'XTickLabel', listROIs );
set(gca,'YTick',1:numROIs,'YTickLabel', listROIs );
title('AUC Old for retrieving pairs of sources')
set(gcf,'position',[100,100,1700,500])
saveas(gcf,['figures' filesep 'pairsAUC'],'png')

figure;
subplot(1,2,1);
imagesc(squeeze(mean(energyNew))); colorbar; caxis([0 1]);
title('Energy New for retrieving pair of sources')
set(gca,'XTick',1:numROIs,'XTickLabel', listROIs );
set(gca,'YTick',1:numROIs,'YTickLabel', listROIs );
subplot(1,2,2);
imagesc(squeeze(mean(energyOld))); colorbar; caxis([0 1]);
set(gca,'XTick',1:numROIs,'XTickLabel', listROIs );
set(gca,'YTick',1:numROIs,'YTickLabel', listROIs );
title('Energy Old for retrieving pairs of sources')    
set(gcf,'position',[100,100,1700,500])
saveas(gcf,['figures' filesep 'pairsEnergy'],'png')

figure;
subplot(1,2,1);
imagesc(squeeze(mean(mseNew))); colorbar; caxis([0 1]);
title('MSE New for retrieving pair of sources')
set(gca,'XTick',1:numROIs,'XTickLabel', listROIs );
set(gca,'YTick',1:numROIs,'YTickLabel', listROIs );
subplot(1,2,2);
imagesc(squeeze(mean(mseOld))); colorbar; caxis([0 1]);
set(gca,'XTick',1:numROIs,'XTickLabel', listROIs );
set(gca,'YTick',1:numROIs,'YTickLabel', listROIs );
title('MSE Old for retrieving pairs of sources')
set(gcf,'position',[100,100,1700,500])
saveas(gcf,['figures' filesep 'pairsMSE'],'png')
