clearvars;close all;
% manipulate the number of sources (picked randomly) per simulation
% ALWAYS unilateral sources 
% only 2 windows: baseline and sources

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
noiseLevel = 10; % 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, reg constant, hyper param in min norm

% nbSbjToInclude =[2 8 20 50];
numSubs = 20;

totBoot = 3;

for totROI=1:numROIs/2
    
    for repBoot=1:totBoot
        fprintf('ROI%d bootstrap %d\n',totROI,repBoot)
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
        if mod(repBoot,2) % only left sources
            uniROI = 1:2:numROIs;
        else % only right sources
            uniROI = 2:2:numROIs;
        end
        ac_sources = randsample(uniROI,totROI);
        source = listROIs(ac_sources);
        
        % amplitude (1 to 10) and time function is different for each
        % source but the same for all sbj for a given bootstrap
        % 1-45 = baseline
        % 46-90 = all sources
        srcERP = zeros(numROIs,45*2);
        srcERP(:,46:90) = createSourceERP(numROIs,ac_sources,[]);
        
        % ERP & baseline timewindow
        timeBase = 1:45;
        winERP = 46:90;
        
        %%% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
        Y_noise = Y;
        Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
        
        for iSub=1:numSubs
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
            for ss=1:length(ac_sources)
                % note that if there is overlapping index (same idx for 2
                % ROIs), the value in sourceData will be of the latest
                % source
                sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
            end
            y_stim = fullFwd{iSub} * sourceData;
            % add noise
            [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
            % to keep the same SNR for the 2 Y, need to compute noise for
            % the 2 Y separately as it is based on the variance of the signal
            Y(iSub,:,:) = y_stim + noisy_data;
        end
        
        %%% Use average reference for centering Y
        for iSub=1:numSubs
            Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
        end
        
        
        %% compute minimum norm
        % min_norm on average data: get beta values for each ROI over time
        [betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
        [betaAverageUni, lambdaUni] = minNormFast_lcurve(avMap(:,uniROI), squeeze(mean(Y_avg,1)));
        
        regionWhole = zeros(numSubs,numROIs,length(srcERP));
        regionROI = zeros(numSubs,numROIs,length(srcERP));
        regionWholeUni = zeros(numSubs,numROIs/2,length(srcERP));
        regionROIUni = zeros(numSubs,numROIs/2,length(srcERP));
        
        for iSub=1:numSubs
            % minimum_norm: on the 20484 indexes per sbj
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
            regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
            
            % minimum_norm: on the 20484 indexes per sbj
            [betaWholeUni,lambdaWholeUni] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            [betaROIUni, lambdaROIUni] = minNormFast([roiFwd{iSub,uniROI}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
     
            % beta values are for the indexes, but I want it per ROI
            % get the number of indexes per ROI for this subj
            rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),uniROI,'uni',false));
            % get the range
            range = [0 cumsum(rangeROI)]; % cumulative sum of elements
            % SUM (not average) the beta values per ROI (=across the indexes)
            regionROIUni(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROIUni(range(x)+1:range(x+1), :)),1:numROIs/2,'uni',false)');
            
            % need to find the indexes for whole brain -> use idxROIfwd
            regionWholeUni(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWholeUni(idxROIfwd{iSub,x},:)),uniROI,'uni',false)');
            
            % feed ROI per sbj instead of mesh
            sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
            [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            sbjROIUni = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),uniROI,'uni',false));
            [betaROIinUni(iSub,:,:), lambdaGridMinNormROIinUni] = minNormFast(sbjROIUni, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole,1));
        retrieveROI = squeeze(mean(regionROI,1));
        retrieveROIin = squeeze(mean(betaROIin,1));
        retrieveWholeUni = squeeze(mean(regionWholeUni,1));
        retrieveROIUni = squeeze(mean(regionROIUni,1));
        retrieveROIinUni = squeeze(mean(betaROIinUni,1));
        
        % save simulation
        simulUni(repBoot,totROI).listROIs = listROIs;
        simulUni(repBoot,totROI).listSub = listSub;
        simulUni(repBoot,totROI).source = source ;
        simulUni(repBoot,totROI).winERP = winERP;
        simulUni(repBoot,totROI).srcERP = srcERP;
        simulUni(repBoot,totROI).srcERPUni = srcERP(uniROI,:);
        simulUni(repBoot,totROI).data = Y_avg;
        simulUni(repBoot,totROI).beta(1,:,:) = betaAverage;
        simulUni(repBoot,totROI).beta(2,:,:) = retrieveWhole;
        simulUni(repBoot,totROI).beta(3,:,:) = retrieveROI;
        simulUni(repBoot,totROI).beta(4,:,:) = retrieveROIin;
        simulUni(repBoot,totROI).betaUni(1,:,:) = betaAverageUni;
        simulUni(repBoot,totROI).betaUni(2,:,:) = retrieveWholeUni;
        simulUni(repBoot,totROI).betaUni(3,:,:) = retrieveROIUni;
        simulUni(repBoot,totROI).betaUni(4,:,:) = retrieveROIinUni;        
        
    end % boot
end % nb of ROI

save('simulOutput/simulRandUni.mat','simulUni')

