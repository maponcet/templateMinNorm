clearvars;close all;
% simulation: V1+MT
% for a given simulation, amplitude and time function is different for each
% source but the same across participants (with different fwd models)
% simulation consistent with retinotopy: L&R sources are the same
% (= assumes full field stimulation)
% simulation on MESH
% minimum_norm: done on a)whole brain, b)only 18 ROIs, c)average50ROIs
% test the effect of SNR for the same nb of sbj only with SSVEP

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.05 0.1 0.2 0.3 0.4 0.5]; % noise level 10% = if the signal is 10 then the noise is 10*10, for 0.5 S=50/100
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm
numCols = 5; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

% nbSbjToInclude =[1 2 5 10 20 30 40 50];
nbSbjToInclude = [20 50];
totBoot = 20; % nb of bootstrap

% initialise variables
ySNR_SSVEP = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
yloSNR_SSVEP = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
aucAveSSVEP = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyAveSSVEP = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseAveNormSSVEP = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));

for totSbj=1:length(nbSbjToInclude)
    numSubs = nbSbjToInclude(totSbj);

for level=1:length(SNRlevel)
        noiseLevel = SNRlevel(level);
        
    for repBoot =1:totBoot
        fprintf('SNR%d bootstrap %d \n',level,repBoot);
        
        % list of random sbj with replacement
        listSub = randi(length(dirList),numSubs,1);
        
        %% LOAD FWD
        fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
        for iSub=1:numSubs
            clear fwdMatrix roiInfo
            % fwd file
            load([dataPath dirList(listSub(iSub)).name])
            fullFwd{iSub} = fwdMatrix;            
%             indexROI = cell2mat(arrayfun(@(x) cellfind({roiInfo.name},listROIs{x}),1:length(listROIs),'uni',false));
            % go through each ROI and save the corresponding fwdMesh values
            % corresponding to the indexes of that ROI
            for rr=1:numROIs
                indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
                roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices); 
                % to get roiFwd for one sbj= [roiFwd{iSub,:}]
                % save the index for each ROI 
                idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
            end
        end
        
        %% Simulate sources (sourceERP)
        % amplitude (1 to 10) and time function is different for each 
        % source but the same for all sbj for a given bootstrap
        [srcAmp, srcSSVEP, srcERP,timeERP] = createSourceROI(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
        % ERP baseline timewindow
        timeBase = setdiff(1:size(srcERP,2),timeERP);
        
        %% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj 
        % (using individual fwd model)
        Y_SSVEP = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
        Y_SSVEPlo = Y_SSVEP;
        Y_noiseSSVEP=Y_SSVEP;
%         Y_SSVEPprev = Y;
        
        for iSub=1:numSubs
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
            sourceDataSSVEP = zeros(size(fullFwd{iSub},2) , length(srcSSVEP));
            sourceNoise = sourceData; % no signal, used to compute SNR
            if length([idxROIfwd{iSub,ac_sources}]) ~= length(unique([idxROIfwd{iSub,ac_sources}]))
                fprintf('Overlapping source indexes S%d \n',listSub(iSub));
            end
            for ss=1:length(ac_sources)
                % note that if there is overlapping index (same idx for 2
                % ROIs), the value in sourceData will be of the latest
                % source
                sourceDataSSVEP(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcSSVEP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
            end
            % multiply fwd (128*20484) with the activated idx over time
            % (sourceData of 20484*90) and obtain Y elec x time
            y_stimSSVEP = fullFwd{iSub} * sourceDataSSVEP;
%             y_noise = fullFwd{iSub} * sourceNoise;
            % add noise
            % to keep the same SNR for the 2 Y, need to compute noise for 
            % the 2 Y separately as it is based on the variance of the signal
            [noisy_dataSSVEP,~] = add_noise_with_SNR( y_stimSSVEP , noiseLevel ); 
            Y_SSVEP(iSub,:,:) = y_stimSSVEP + noisy_dataSSVEP;
%             Y_SSVEPprev(iSub,:,:) = y_stimSSVEP + noisy_data;
%             Y_noise(iSub,:,:) = noisy_data;
            Y_noiseSSVEP(iSub,:,:) = noisy_dataSSVEP;
%             Ypure(iSub,:,:) = y_stim;
%             Y2(iSub,:,:) = y_stim + noisy_dataSSVEP;
        end
        % check SNR
        currSNR = zeros(1,numSubs);
        for iSub=1:numSubs
            currSNR(iSub) = (rms(rms(Y_SSVEP(iSub,:,:)))/rms(rms(Y_noiseSSVEP(iSub,:,:)))) ^2 -1;
        end       
        ySNR_SSVEP(repBoot,totSbj,level) = mean(currSNR);    
        
        %         figure;plotContourOnScalp(squeeze(Y(1,:,45)),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
        
        %% center Y, fullFwd, roiFwd across electrodes?? Sbj??
        % use average reference for centering the data
        for iSub=1:numSubs
%             Y2(iSub,:,:) = bsxfun(@minus,squeeze(Y2(iSub,:,:)), mean(squeeze(Y2(iSub,:,:))));
            Y_SSVEP(iSub,:,:) = bsxfun(@minus,squeeze(Y_SSVEP(iSub,:,:)), mean(squeeze(Y_SSVEP(iSub,:,:))));
        end
        
        %%  reduce dimension data (Ylo)
        % =PCA denoised version of Y (denoised by truncation of the SVD)
        %%%%%% per sbj or alltogether?
        for iSub=1:numSubs
%             [u1, s1, v1] = svd(squeeze(Y2(iSub,:,:)));
%             Y2lo(iSub,:,:) = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
            [u2, s2, v2] = svd(squeeze(Y_SSVEP(iSub,:,:)));
            Y_SSVEPlo(iSub,:,:) = u2(:,1:numCols)*s2(1:numCols,1:numCols)*v2(:, 1:numCols)';
        end
        % compute SNR after reduction... without -1 because i only have
        % signal now???
        yloSNR = zeros(1,numSubs);
        for iSub=1:numSubs
            yloSNR(iSub) = (rms(rms(Y_SSVEPlo(iSub,:,:)))/rms(rms(Y_noiseSSVEP(iSub,:,:)))) ^2 ;          
        end
        yloSNR_SSVEP(repBoot,totSbj,level) = mean(yloSNR);    
        
%        % stack sbj vertically
%        tmpY = permute(Y,[2 1 3]);
%        tmpY2 = reshape(tmpY,[size(fullFwd{1},1)*numSubs,length(srcERP)]);
%        [u1, s1, v1] = svd(tmpY2);
%        YloStack = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
%        tmpYs = permute(Y_SSVEP,[2 1 3]);
%        tmpY2s = reshape(tmpYs,[size(fullFwd{1},1)*numSubs,length(srcERP)]);
%        [u1, s1, v1] = svd(tmpY2s);
%        Y_SSVEPloStack = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
       
        %% compute minimum norm
%         regionWhole = zeros(numSubs,numROIs,length(srcERP));
%         regionROI = zeros(numSubs,numROIs,length(srcERP));        
        regionWhole_SSVEP = zeros(numSubs,numROIs,length(srcERP));
        regionROI_SSVEP = zeros(numSubs,numROIs,length(srcERP));
        % min_norm on average data: get beta values for each ROI over time 
        [betaAverageSSVEP, ~, lambdaAverageSSVEP] = minimum_norm(avMap, squeeze(mean(Y_SSVEPlo,1)), nLambdaRidge);
%         unstackedY = reshape(YloStack,size(fullFwd{1},1),numSubs,size(YloStack,2));
%         [betaAverageStack, ~, lambdaAverageStack] = minimum_norm(avMap, squeeze(mean(unstackedY,2)), nLambdaRidge);
%         unstackedSSVEP = reshape(Y_SSVEPloStack,size(fullFwd{1},1),numSubs,size(Y_SSVEPloStack,2));
%         [betaAverageSSVEP, ~, lambdaAverageSSVEP] = minimum_norm(avMap, squeeze(mean(unstackedSSVEP,2)), nLambdaRidge);
%         for iSub=1:numSubs
%             % regular minimum_norm: on the 20484 indexes per sbj
% %             [betaWhole, ~, lambdaWhole] = minimum_norm(fullFwd{iSub}, squeeze(Ylo(iSub,:,:)), nLambdaRidge);
%             [betaWhole_SSVEP, ~, lambdaWhole_SSVEP] = minimum_norm(fullFwd{iSub}, squeeze(Y_SSVEPlo(iSub,:,:)), nLambdaRidge);
% %             [betaWhole, ~, lambdaWhole] = minimum_norm(fullFwd{iSub}, squeeze(unstackedY(:,iSub,:)), nLambdaRidge);
%             % min_norm on only the ROI indexes per sbj
% %             [betaROI, ~, lambdaROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Ylo(iSub,:,:)), nLambdaRidge);
%             [betaROI_SSVEP, ~, lambdaROI_SSVEP] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_SSVEPlo(iSub,:,:)), nLambdaRidge);
% %             [betaROIStack, ~, lambdaROIStack] = minimum_norm([roiFwd{iSub,:}], squeeze(unstackedY(:,iSub,:)), nLambdaRidge);
%             
%             % beta values are for the indexes, but I want it per ROI
%             % get the number of indexes per ROI for this subj
%             rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
%             % get the range
%             range = [0 cumsum(rangeROI)]; % cumulative sum of elements
%             % average the beta values per ROI (=across the indexes)
% %             regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); 
%             regionROI_SSVEP(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROI_SSVEP(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); 
% %             regionROI_Stack(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROIStack(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); 
%             % need to find the indexes for whole brain
%             regionWhole_SSVEP(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole_SSVEP(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%             % create brain resp for each sbj
%             % YhatWhole(iSub,:,:) = fullFwd{iSub}*betaWhole(iSub,:,:);
%         end
        % average across subj
%         retrieveWhole = squeeze(mean(regionWhole,1));
%         retrieveROI = squeeze(mean(regionROI,1));
%         retrieveROI_SSVEP = squeeze(mean(regionROI_SSVEP,1));
%         retrieveWhole_SSVEP = squeeze(mean(regionWhole_SSVEP,1));
%         retrieveROI_Stack = squeeze(mean(regionROI_Stack,1));
        
        %% compute auc, mse, relative energy using average signal in rois
        % do for all the min norm outputs
%         [aucWhole(repBoot,totSbj), energyWhole(repBoot,totSbj),mseWhole(repBoot,totSbj),...
%             mseWholeNorm(repBoot,totSbj)] = computeMetrics(retrieveWhole,ac_sources,srcERP);
%         [aucROI(repBoot,totSbj), energyROI(repBoot,totSbj),mseROI(repBoot,totSbj),...
%             mseROINorm(repBoot,totSbj)] = computeMetrics(retrieveROI,ac_sources,srcERP);
        [aucAveSSVEP(repBoot,totSbj,level), energyAveSSVEP(repBoot,totSbj,level),...
            mseAveNormSSVEP(repBoot,totSbj,level)] = computeMetrics(betaAverageSSVEP,ac_sources,srcSSVEP);
%         [aucAveStack(repBoot,totSbj), energyAveStack(repBoot,totSbj),...
%             mseAveStack(repBoot,totSbj)] = computeMetrics(betaAverageStack,ac_sources,srcERP);
%         [aucROI_SSVEP(repBoot,level), energyROI_SSVEP(repBoot,level),...
%             mseROINorm_SSVEP(repBoot,level)] = computeMetrics(retrieveROI_SSVEP,ac_sources,srcSSVEP);
%         [aucWhole_SSVEP(repBoot,level), energyWhole_SSVEP(repBoot,level),...
%             mseWholeNorm_SSVEP(repBoot,level)] = computeMetrics(retrieveWhole_SSVEP,ac_sources,srcSSVEP);
%         [aucROI_Stack(repBoot,totSbj), energyROI_Stack(repBoot,totSbj),...
%             mseROI_Stack(repBoot,totSbj)] = computeMetrics(retrieveROI_Stack,ac_sources,srcERP);        
    end % end boostrap

    
end % go through nb of SNR level to include

end % nbSbjToInclude

figure;
subplot(2,2,1);hold on;hold on
for ss=1:2
errorbar(SNRlevel,squeeze(mean(aucAveSSVEP(:,ss,:))),squeeze(std(aucAveSSVEP(:,ss,:))))
end
xlabel('SNR')
ylabel('AUC')
legend('N20','N50')
subplot(2,2,2);hold on;hold on
for ss=1:2
errorbar(SNRlevel,squeeze(mean(energyAveSSVEP(:,ss,:))),squeeze(std(energyAveSSVEP(:,ss,:))))
end
xlabel('SNR')
ylabel('energy')
legend('N20','N50')
subplot(2,2,3);hold on;hold on
for ss=1:2
errorbar(SNRlevel,squeeze(mean(mseAveNormSSVEP(:,ss,:))),squeeze(std(mseAveNormSSVEP(:,ss,:))))
end
xlabel('SNR')
ylabel('normMSE')
legend('N20','N50')
subplot(2,2,4);hold on
% for ss=1:2
    errorbar(SNRlevel,squeeze(mean(ySNR_SSVEP(:,ss,:))),squeeze(std(ySNR_SSVEP(:,ss,:))))
    errorbar(SNRlevel,squeeze(mean(yloSNR_SSVEP(:,ss,:))),squeeze(std(yloSNR_SSVEP(:,ss,:))))
% end
xlabel('SNR')
ylabel('SNR')
legend('beforeICA','afterICA')
saveas(gcf,['figures' filesep 'SNRtest-SSVEP'],'png')
