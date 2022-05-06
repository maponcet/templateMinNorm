clearvars;close all;
% simulate ERP from outside the 18 ROIs by picking some sources that are
% >15000

% retrieve the sources using different methods (template, ROI, whole brain)
% Each simulation uses a different set of sbj and ERP but the same one is
% tested across different levels of noise and number of sbj

addpath(([pwd filesep 'subfunctions']))
addpath('/Users/marleneponcet/Documents/Git/ssvepTesting/svndlCopy/') % for plotEGI
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = 200; % [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm

nbSbjToInclude = 20; % [2 8 20 50];

totBoot = 10;

for repBoot=1:totBoot
   
    for totSbj=1:length(nbSbjToInclude) 
        numSubs = nbSbjToInclude(totSbj);
        fprintf('N%d bootstrap %d\n',numSubs,repBoot)
        
        % list of random sbj with replacement
        listSub = randi(length(dirList),numSubs,1);
        
        %% LOAD FWD
        fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
        for iSub=1:numSubs
            clear fwdMatrix roiInfo
            % fwd file
            load([dataPath dirList(listSub(iSub)).name])
            fullFwd{iSub} = fwdMatrix;
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
        
        for level=1:length(SNRlevel)
            noiseLevel = SNRlevel(level);
            
            %%% Simulate scalp activity (Y)
            % define time windows
            timeBase = 1:5;winERP = 6:10; fullWin = [timeBase winERP];
            
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{1},2) , length(fullWin));
            Y = zeros(numSubs,size(fullFwd{1},1),length(fullWin));
            Y_avg = zeros(numSubs,size(fullFwd{1},1),length(fullWin));
 
            % pick 100 random active sources all >15000
            % if too many sources or random for each sbj then just get same
            % average topo each simul
            activeROI = randi([15000 size(fwdMatrix,2)],100,1);
            sourceData(activeROI,winERP) = 10;
            for iSub=1:numSubs
                y_stim = fullFwd{iSub} * sourceData;
                % add noise
                [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
                Y(iSub,:,:) = y_stim + noisy_data;
            end

            %%% Use average reference for centering Y
            for iSub=1:numSubs
                Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
            end

 
            %% compute minimum norm
%             regionWhole = zeros(numSubs,numROIs,length(srcERP));
%             regionROI = regionWhole;
%             regionROILC = regionWhole;
%             regionWholeLC = regionWhole;
%             betaROIin = regionWhole;
%             betaROIinLC = regionWhole;
            
            % min_norm on average data: get beta values for each ROI over time
            [betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
            
            for iSub=1:numSubs
                % regular minimum_norm: on the 20484 indexes per sbj
                [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
                [betaWholeLC, lambdaWholeLC] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
                
                [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
                [betaROILC, lambdaROILC] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
                
                % beta values are for the indexes, but needs it per ROI
                % get the number of indexes per ROI for this subj
                rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
                % get the range
                range = [0 cumsum(rangeROI)]; % cumulative sum of elements
                % SUM (not average) the beta values per ROI (=across the indexes)
                regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
                regionROILC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
                
                % need to find the indexes for whole brain -> use idxROIfwd
                % (no need to get the range)
                regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
                regionWholeLC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
                
            end
            % average across subj
            retrieveWhole = squeeze(mean(regionWhole,1));
            retrieveROI = squeeze(mean(regionROI,1));
            retrieveWholeLC = squeeze(mean(regionWholeLC,1));
            retrieveROILC = squeeze(mean(regionROILC,1));
            
            figure;subplot(3,2,1);imagesc(betaAverage/max(betaAverage(:)));caxis([0 1]);title('template')
            subplot(3,2,3);imagesc(retrieveWhole/max(retrieveWhole(:)));caxis([0 1]);title('whole')
            subplot(3,2,4);imagesc(retrieveWholeLC/max(retrieveWholeLC(:)));caxis([0 1])
            subplot(3,2,5);imagesc(retrieveROI/max(retrieveROI(:)));caxis([0 1]);title('ROI')
            subplot(3,2,6);imagesc(retrieveROILC/max(retrieveROILC(:)));caxis([0 1])
            saveas(gcf,['figures/outROIsim' num2str(repBoot) 'betas'],'png')
            
            %%% plot topographies at t=8
            predTemp = avMap * betaAverage;
            predWhole = avMap * retrieveWhole;
            predWholeLC = avMap * retrieveWholeLC;
            predROI = avMap * retrieveROI;
            predROILC = avMap * retrieveROILC;
            
            figure;set(gcf, 'Position', [0 0 900 700])
            subplot(3,2,1);plotOnEgi(squeeze(mean(Y_avg(:,:,8),1))');title('data')
            subplot(3,2,2);plotOnEgi(predTemp(:,8));title('template')
            subplot(3,2,3);plotOnEgi(predWhole(:,8));title('whole')
            subplot(3,2,4);plotOnEgi(predWholeLC(:,8));title('wholeLC')
            subplot(3,2,5);plotOnEgi(predROI(:,8));title('ROI')
            subplot(3,2,6);plotOnEgi(predROILC(:,8));title('ROILC')
            saveas(gcf,['figures/outROIsim' num2str(repBoot)],'png')
        end % noise
        
    end % sbj
end % boot

