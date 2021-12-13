clearvars;close all;
% simulate one area and compute amount of leakage (crosstalk) with other
% areas

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
crossTalkTemplate = zeros(totBoot,numROIs,numROIs);
crossTalkMinNorm = zeros(totBoot,numROIs,numROIs);
crossTalkTemplate_ERP = zeros(totBoot,numROIs,numROIs);
crossTalkMinNorm_ERP = zeros(totBoot,numROIs,numROIs);

for seedRoi = 1:length(listROIs)
    
%     % simulated signal
%     signalROIs = listROIs(seedRoi);
%     % find the ROI index which corresponds to the one in signalROIs
%     ac_sources = cell2mat(arrayfun(@(x) cellfind(avMap.roiNames,signalROIs{x}),1:length(signalROIs),'uni',false));
    
    for repBoot =1:totBoot
        fprintf('sim%d bootstrap %d \n',seedRoi,repBoot);
        clear Y source signal Ylo betaReg betaMinNorm
        
        % list of random sbj with replacement
        listSub = randi(length(dirList),numSubs,1);
        
        %% GET ROIs
        % do for each sbj separately but on 18 ROIs, not mesh
        % so first need to get the ROI per sbj
        roiMap=zeros(128,numROIs,numSubs);
        roiStacked = [];
        for iSubj=1:numSubs
            clear fwdMatrix roiInfo
            % fwd file
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
        % amplitude and time function is different for each source but
        % the same across participants for a given bootstrap
        x = 0 : pi / 45 : 2 * pi-pi/45;
        sourceValOverTime = zeros(numROIs, length(x) );
        sourceValOverTimeERP = zeros(numROIs, length(x) );
        % random amplitude between 1 and 10
        sourceAmplitude = 1 + randi(90)/10;
        % create SSVEP waveform with random phase delay and random
        % coef for each source
        y = rand * cos( 4 * x - rand * pi ) + ...
            rand * cos( 8 * x - rand * pi ) + ...
            rand * cos( 12 * x - rand * pi ) + ...
            rand * cos( 16 * x - rand * pi );
        % create ERP waveform (use cos envelop instead?)
        yERP = y.*gaussmf(x,[pi/5 pi]);
        % multiply waveform by the source amplitude (ROIxtime)
        sourceValOverTime( seedRoi , : ) = y * sourceAmplitude;
        sourceValOverTimeERP( seedRoi , : ) = yERP * sourceAmplitude;
            
        % Simul EEG signal from sourceValOverTime for each sbj (different
        % ROI map)
        for iSub=1:numSubs
            y_stim = roiMap(:,:,iSub) * sourceValOverTime;
            y_stimERP = roiMap(:,:,iSub) * sourceValOverTimeERP;
            % add noise
            [noisy_data,noise_var] = add_noise_with_SNR( y_stim , SNRlevel );
            Y(128*(iSub-1)+1:128*iSub,:) = y_stim + noisy_data;
            [noisy_data,noise_var] = add_noise_with_SNR( y_stimERP , SNRlevel );
            Yerp(128*(iSub-1)+1:128*iSub,:) = y_stimERP + noisy_data;
        end
        
        %%  reduce dimension data (Ylo)
        % =PCA denoised version of Y (denoised by truncation of the SVD)
        n = numel(Y);
        [u1, s1, v1] = svd(Y);
        Ylo = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
        % compute average EEG
        unstackedData = reshape(Ylo,128,numSubs,size(Ylo,2));
        grandMeanData = squeeze(mean(unstackedData,2));
        % do the same for ERP simul
        [u1ERP, s1ERP, v1ERP] = svd(Yerp);
        YloERP = u1ERP(:,1:numCols)*s1ERP(1:numCols,1:numCols)*v1ERP(:, 1:numCols)';
        unstackedDataERP = reshape(YloERP,128,numSubs,size(YloERP,2));
        grandMeanDataERP = squeeze(mean(unstackedDataERP,2));
            
        %% compute minimum norm
        [betaReg, ~, lambdaReg] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
        [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(roiStacked, Ylo, nLambdaRidge);
        [betaRegERP, ~, lambdaRegERP] = minimum_norm(avMap.activity, grandMeanDataERP, nLambdaRidge);
        [betaMinNormERP, ~, lambdaMinNormERP] = minimum_norm(roiStacked, YloERP, nLambdaRidge);
        % betaMinNorm is numSubs*numROIs x time so need to extract betaMinNorm of a single sbj
        % for each ROI separately then average across sbj
        regionTmp = arrayfun(@(x) betaMinNorm(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
        regionMinNorm =reshape(cell2mat(regionTmp),numROIs,size(Ylo,2),numSubs);
        regionMinNormAvg = mean(regionMinNorm,3);
        regionTmpERP = arrayfun(@(x) betaMinNormERP(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
        regionMinNormERP =reshape(cell2mat(regionTmpERP),numROIs,size(Ylo,2),numSubs);
        regionMinNormAvgERP = mean(regionMinNormERP,3);
        
        %% amount of crosstalk
        % leakage: amplitude signal in all areas, the "true" one is used
        % for normalisation (=1). Use RMS
        normTerm = rms(betaReg(seedRoi,:));
        normTermMinNorm = rms(regionMinNormAvg(seedRoi,:));
        normTermERP = rms(betaRegERP(seedRoi,:));
        normTermMinNormERP = rms(regionMinNormAvgERP(seedRoi,:));
        for iRoi = 1:numROIs
            crossTalkTemplate(repBoot,seedRoi,iRoi) = rms(betaReg(iRoi,:)) / normTerm;
            crossTalkMinNorm(repBoot,seedRoi,iRoi) = rms(regionMinNormAvg(iRoi,:)) / normTermMinNorm;
            crossTalkTemplate_ERP(repBoot,seedRoi,iRoi) = rms(betaRegERP(iRoi,:)) / normTermERP;
            crossTalkMinNorm_ERP(repBoot,seedRoi,iRoi) = rms(regionMinNormAvgERP(iRoi,:)) / normTermMinNormERP;
        end
        
    end % end boostrap
    
end % different activated sources

figure;
subplot(1,2,1);imagesc(squeeze(mean(crossTalkTemplate)));colorbar;caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs)   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs)  
ylabel('seedArea');xlabel('predictArea')
title('templateBased')
subplot(1,2,2);imagesc(squeeze(mean(crossTalkMinNorm)));colorbar;caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs)   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs)   
ylabel('seedArea');xlabel('predictArea')
title('Standard MinNorm')
set(gcf,'position',[100,100,1800,800])
saveas(gcf,['figures' filesep 'crossTalk'],'png')

figure;
subplot(1,2,1);imagesc(squeeze(mean(crossTalkTemplate_ERP)));colorbar;caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs)   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs)  
ylabel('seedArea');xlabel('predictArea')
title('templateBasedERP')
subplot(1,2,2);imagesc(squeeze(mean(crossTalkMinNorm_ERP)));colorbar;caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs)   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs)   
ylabel('seedArea');xlabel('predictArea')
title('Standard MinNormERP')
set(gcf,'position',[100,100,1800,800])
saveas(gcf,['figures' filesep 'crossTalkERP'],'png')
