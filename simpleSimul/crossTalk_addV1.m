clearvars;close all;
% simulate one area and compute amount of leakage (crosstalk) with other
% areas

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
numSubs = length(dirList);

% some parameters
nLambdaRidge = 20; % for calculating minimum_norm, reg constant, hyper param in min norm
SNRlevel = [1 10 200]; % noise level
totBoot = 30; % nb of bootstrap
crossTalkTemplate = zeros(totBoot,numROIs,numROIs);
crossTalkWhole = zeros(totBoot,numROIs,numROIs);
crossTalkROI = zeros(totBoot,numROIs,numROIs);
crossTalkROIin = crossTalkTemplate;
crossTalkTemplateBest= crossTalkTemplate;

crossTalkNormTemplate = crossTalkTemplate;
crossTalkNormWhole = crossTalkTemplate;
crossTalkNormROI = crossTalkTemplate;
crossTalkNormROIin = crossTalkTemplate;
crossTalkNormTemplateBest = crossTalkTemplate;


%% LOAD FWD - just keep the same 50 participants (ie no randomisation)
fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
for iSub=1:numSubs
    clear fwdMatrix roiInfo
    % fwd file
    load([dataPath dirList(iSub).name])
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
        
for nn=1:length(SNRlevel)  
    noise = SNRlevel(nn);
for seedRoi = 2:numROIs
    
    activeROIs = [{'V1-L'}, listROIs(seedRoi)];
    ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));
    
    for repBoot =1:totBoot
        fprintf('seed%d bootstrap %d noise %d\n',seedRoi,repBoot,nn);
        clear Y source signal Ylo betaReg betaMinNorm

        %% Simulate sources 
        srcERP = zeros(numROIs, 20 );
        srcERP(ac_sources,:) = 1;        
        
        %% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj 
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),size(srcERP,2));   
        Y_avg = Y;
        for iSub=1:numSubs
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{iSub},2) , size(srcERP,2));
            sourceNoise = sourceData; % no signal, used to compute SNR
            for ss=1:length(ac_sources)
                % note that if there is overlapping index (same idx for 2
                % ROIs), the value in sourceData will be of the latest
                % source
                sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
            end
            % multiply fwd (128*20484) with the activated idx over time
            % (sourceData of 20484*90) and obtain Y elec x time
            y_stim = fullFwd{iSub} * sourceData;
            % add gaussian noise on electrodes
            sig = rms(y_stim(:)); % get amount of signal
            noisy_data = sig/sqrt(noise) * randn(size(y_stim,1),size(y_stim,2));
            Y(iSub,:,:) = y_stim + noisy_data;
%             % check SNR
%             allSig = Y(iSub,:,:);
%             (rms(allSig(:))/rms(noisy_data(:)))^2 - 1
        end
    
    
        %% Use average reference for centering Y
        %%% that is: substract the average electrode activity at each time point
        % this is done by bsxfun which applies element-wise substraction (the 90
        % averages across electrodes)
        for iSub=1:numSubs
            Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
        end

        %% compute minimum norm    
        % min_norm: get beta values for each ROI over time 
        [betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
%         [betaLCFUN, betaCurv, betaBest, lambdaLCFUN, lambdaCurv, lambdaBest, ...
%             lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(Y_avg,1)),srcERP);
        
%         regionWhole = zeros(numSubs,numROIs,size(srcERP,2));
%         regionROI = zeros(numSubs,numROIs,size(srcERP,2));
%         betaROIin = regionWhole;
%         betaROIinLC = regionWhole;
%         
%         for iSub=1:numSubs
%             % regular minimum_norm: on the 20484 indexes per sbj
%             [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%             [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%      
%             % beta values are for the indexes, but I want it per ROI
%             % get the number of indexes per ROI for this subj
%             rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
%             % get the range
%             range = [0 cumsum(rangeROI)]; % cumulative sum of elements
%             % SUM (not average) the beta values per ROI (=across the indexes)
%             regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%             
%             % need to find the indexes for whole brain -> use idxROIfwd
%             % (no need to get the range)
%             regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%             
%             % feed ROI per sbj instead of mesh
%             sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
%             [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%         end
%         % average across subj
%         retrieveWhole = squeeze(mean(regionWhole,1));
%         retrieveROI = squeeze(mean(regionROI,1));
%         retrieveROIin = squeeze(mean(betaROIin,1));

        %% amount of crosstalk
        % leakage: amplitude signal in all areas (NOT normalised)
        crossTalkTemplate(repBoot,seedRoi,:) = rms(betaAverage,2) ;
%         crossTalkWhole(repBoot,seedRoi,:) = rms(retrieveWhole,2) ;
%         crossTalkROI(repBoot,seedRoi,:) = rms(retrieveROI,2) ;
%         crossTalkROIin(repBoot,seedRoi,:) = rms(retrieveROIin,2) ;
%         crossTalkTemplateBest(repBoot,seedRoi,:) = rms(betaBest,2) ;   
        % normalise
        crossTalkNormTemplate(repBoot,seedRoi,:) = rms(betaAverage,2) / max(rms(betaAverage,2)) ;
%         crossTalkNormWhole(repBoot,seedRoi,:) = rms(retrieveWhole,2) / max(rms(retrieveWhole,2));
%         crossTalkNormROI(repBoot,seedRoi,:) = rms(retrieveROI,2) / max(rms(retrieveROI,2));
%         crossTalkNormROIin(repBoot,seedRoi,:) = rms(retrieveROIin,2) / max(rms(retrieveROIin,2));
%         crossTalkNormTemplateBest(repBoot,seedRoi,:) = rms(betaBest,2) / max(rms(betaBest,2));          
        
    end % end boostrap
    
end % different activated sources

save(['simulOutput/crossTalkV1N' num2str(SNRlevel(nn)) '.mat'],'listROIs',...
    'crossTalkTemplate','crossTalkNormTemplate');

% save(['simulOutput/crossTalkN' num2str(SNRlevel(nn)) '.mat'],'listROIs',...
%     'crossTalkTemplate','crossTalkWhole','crossTalkROI','crossTalkROIin','crossTalkTemplateBest',...
%     'crossTalkNormTemplate','crossTalkNormWhole','crossTalkNormROI','crossTalkNormROIin','crossTalkNormTemplateBest');
end        
   

