clearvars;close all;

addpath(genpath([pwd filesep 'subfunctions']))

% path for the different fwd models
mm(1).dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI32/';
mm(1).dirList = dir([mm(1).dataPath 'forward*']);
mm(2).dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI64/';
mm(2).dirList = dir([mm(2).dataPath 'forward*']);
mm(3).dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/';
mm(3).dirList = dir([mm(3).dataPath 'forward*']);
mm(4).dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI256/';
mm(4).dirList = dir([mm(4).dataPath 'forward*']);
nbModels = length(mm);

dirModel = '/Users/marleneponcet/Documents/Git/templateMinNorm/createTemplate/templates/';
load([dirModel 'averageMapEGI32.mat']); 
mm(1).avMap = avMap;
load([dirModel 'averageMapEGI64.mat']); 
mm(2).avMap = avMap;
load([dirModel 'averageMapEGI128.mat']); % load average map of ROIs (128 elec x 18 ROIs)
mm(3).avMap = avMap;
load([dirModel 'averageMapEGI256.mat']); 
mm(4).avMap = avMap;

numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
% SNRlevel = [0.1  10  10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

numSubs = 50;

totBoot = 30;


for repBoot=1:totBoot
    
    % list of random sbj with replacement
    listSub = randi(numSubs,numSubs,1);

    
    %% Simulate sources
    % amplitude (1 to 10) and time function is different for each
    % source but the same for all sbj for a given bootstrap
    % 1-45 = baseline
    % 46-90 = V1
    % 91-135 = MT
    % 136-180 = MT+V1
    srcERP = zeros(numROIs,45*4); % 45*4 timepoints
    srcERP(:,46:90) = createSourceERP(numROIs,ac_sources(1),ac_sources(3));
    srcERP(:,91:135) = createSourceERP(numROIs,ac_sources(2),ac_sources(4));
    srcERP(:,136:180) = createSourceERP(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
    
    % ERP & baseline timewindow
    timeBase = 1:45;
    winERP = 46:180;
    
    for level=1:length(SNRlevel)
    noiseLevel = SNRlevel(level);

    for sysNb=1:nbModels
        fprintf('bootstrap %d noise %d model %d \n',repBoot,level,sysNb)
        
        %% LOAD FWD
        fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
        for iSub=1:numSubs
            clear fwdMatrix roiInfo
            % fwd file
            load([mm(sysNb).dataPath mm(sysNb).dirList(listSub(iSub)).name])
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
    
        
        %%% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
        Y_noise = Y;
        Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
        
        for iSub=1:numSubs
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
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
            % add noise
            [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
            Y(iSub,:,:) = y_stim + noisy_data;
            Y_noise(iSub,:,:) = noisy_data;
        end
        
        %%% Use average reference for centering Y
        %%% that is: substract the average electrode activity at each time point
        % this is done by bsxfun which applies element-wise substraction (the 90
        % averages across electrodes) - Useless
        for iSub=1:numSubs
            Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
        end

        
        %% compute minimum norm
        [betaAverage, lambda] = minNormFast_lcurve(mm(sysNb).avMap, squeeze(mean(Y_avg,1)));
        
        regionWhole = zeros(numSubs,numROIs,length(srcERP));
        regionROI = regionWhole;
        betaROIin = regionWhole;
        for iSub=1:numSubs
            % regular minimum_norm: on the 20484 indexes per sbj
            [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            
            [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            
            % beta values are for the indexes, but needs it per ROI
            % get the number of indexes per ROI for this subj
            rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
            % get the range
            range = [0 cumsum(rangeROI)]; % cumulative sum of elements
            % SUM (not average) the beta values per ROI (=across the indexes)
            regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
            
            % need to find the indexes for whole brain -> use idxROIfwd
            % (no need to get the range)
            regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
            
            % feed ROI per sbj instead of mesh = "oracle" (best possible recovery)
            sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
            [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole,1));
        retrieveROI = squeeze(mean(regionROI,1));
        retrieveROIin = squeeze(mean(betaROIin,1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % save simulation
        simulSys(repBoot,level,sysNb).listROIs = listROIs;
        simulSys(repBoot,level,sysNb).listSub = listSub;
        simulSys(repBoot,level,sysNb).winERP = winERP;
        simulSys(repBoot,level,sysNb).srcERP = srcERP;
        simulSys(repBoot,level,sysNb).data = Y_avg;
        simulSys(repBoot,level,sysNb).noise = SNRlevel(level);
        simulSys(repBoot,level,sysNb).beta(1,:,:) = betaAverage;
        simulSys(repBoot,level,sysNb).beta(2,:,:) = retrieveWhole;
        simulSys(repBoot,level,sysNb).beta(3,:,:) = retrieveROI;
        simulSys(repBoot,level,sysNb).beta(4,:,:) = retrieveROIin;
        
    end
    end
    
end

save('simulOutput/simulSysV1MT.mat','simulSys','-v7.3')


