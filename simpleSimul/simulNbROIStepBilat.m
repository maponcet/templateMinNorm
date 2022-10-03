clearvars;close all;
% manipulate the number of sources (picked randomly) per simulation
% ALWAYS bilateral sources (if V1 is active, both left & right are)
% only 2 windows: baseline and sources (10 timepoints each)
% source is 1 (step, not ERP; set at 1, not random)

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
noiseLevel = [1 10 200]; 
nLambdaRidge = 20; % for calculating minimum_norm, reg constant, hyper param in min norm

% nbSbjToInclude =[2 8 20 50];
numSubs = 50;

totBoot = 30;

for snr=1:length(noiseLevel)
    noise = noiseLevel(snr);
    
for totROI=1:numROIs/2
    
    for repBoot=1:totBoot
        fprintf('ROI%d bootstrap %d SNR %d\n',totROI,repBoot,noise)
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
        
        %% Simulate sources
        roiActive = randsample(1:2:numROIs,totROI);
        sourceL = listROIs(roiActive);
        sourceR = listROIs(roiActive+1);
        activeROIs = [sourceL,sourceR];
        % find the ROI index corresponding to the activeROIs
        ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));
        
        srcERP = zeros(numROIs,20);
        srcERP(ac_sources,:) = 1;    
        
        %%% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),size(srcERP,2));
        Y_noise = Y;
        Y_avg = zeros(numSubs,size(fullFwd{1},1),size(srcERP,2));
        
        for iSub=1:numSubs
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{iSub},2) , size(srcERP,2));
            for ss=1:length(ac_sources)
                % note that if there is overlapping index (same idx for 2
                % ROIs), the value in sourceData will be of the latest
                % source
                sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
            end
            y_stim = fullFwd{iSub} * sourceData;
            % add gaussian noise on electrodes
            sig = rms(y_stim(:)); % get amount of signal
            noisy_data = sig/sqrt(noise) * randn(size(y_stim,1),size(y_stim,2));
            Y(iSub,:,:) = y_stim + noisy_data;

%             % check SNR
%             allSig = Y(iSub,:,:);
%             (rms(allSig(:))/rms(noisy_data(:)))^2 - 1
        end

        %%% Use average reference for centering Y
        for iSub=1:numSubs
            Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
        end
        
        
        %% compute minimum norm
        % min_norm on average data: get beta values
        % corner computed using either lcfun or curvature (which is chosen
        % here: betaAverage). betabest is with best lambda
        [betaLcFun, betaAverage, betaBest, lambdaLcFun, lambdaCurv, lambdaBest, ...
            lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(Y_avg,1)), srcERP);
        
        regionWhole = zeros(numSubs,numROIs,length(srcERP));
        regionROI = zeros(numSubs,numROIs,length(srcERP));
        betaROIin = regionWhole;
        betaROIinLC = regionWhole;
        
        for iSub=1:numSubs
            % regular minimum_norm: on the 20484 indexes per sbj
            [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
     
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
            
            % feed ROI per sbj instead of mesh
            sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
            [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            [betaROIinLC(iSub,:,:), lambdaGridMinNormROIinLC] = minNormFast_lcurve(sbjROI, squeeze(Y_avg(iSub,:,:)));
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole,1));
        retrieveROI = squeeze(mean(regionROI,1));
        retrieveROIin = squeeze(mean(betaROIin,1));
        retrieveROIinLC = squeeze(mean(betaROIinLC,1));

        % save simulation
        simulBilat(repBoot,totROI).listROIs = listROIs;
        simulBilat(repBoot,totROI).listSub = listSub;
        simulBilat(repBoot,totROI).activeROIs = activeROIs ;
        simulBilat(repBoot,totROI).activeSources = ac_sources ;
        simulBilat(repBoot,totROI).srcERP = srcERP;
        simulBilat(repBoot,totROI).data = Y_avg;
        simulBilat(repBoot,totROI).beta(1,:,:) = betaAverage;
        simulBilat(repBoot,totROI).beta(2,:,:) = retrieveWhole;
        simulBilat(repBoot,totROI).beta(3,:,:) = retrieveROI;
        simulBilat(repBoot,totROI).beta(4,:,:) = retrieveROIin;
        simulBilat(repBoot,totROI).beta(5,:,:) = retrieveROIinLC;
        simulBilat(repBoot,totROI).beta(6,:,:) = betaBest;
        
    end % boot
end % nb of ROI

save(['simulOutput' filesep 'simulStepBilatSNR' num2str(noise) '.mat'],'simulBilat')

end
