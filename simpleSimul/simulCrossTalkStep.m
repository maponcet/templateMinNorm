clearvars;close all;
% simulate one area and compute amount of leakage (crosstalk) with other
% areas

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
numSubs = length(dirList);

% some parameters
nLambdaRidge = 20; % for calculating minimum_norm, reg constant, hyper param in min norm
SNRlevel = 200; % noise level, 200 times more signal than noise
totBoot = 30; % nb of bootstrap
crossTalkTemplate = zeros(totBoot,numROIs,numROIs);
crossTalkWhole = zeros(totBoot,numROIs,numROIs);
crossTalkROI = zeros(totBoot,numROIs,numROIs);
crossTalkROIin = crossTalkTemplate;
crossTalkTemplateBest= crossTalkTemplate;

for seedRoi = 1:numROIs
    
    activeROIs = listROIs(seedRoi);
    ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));
    
    for repBoot =1:totBoot
        fprintf('seed%d bootstrap %d \n',seedRoi,repBoot);
        clear Y source signal Ylo betaReg betaMinNorm
        
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
        
        
        %% Simulate sources ON/OFF
        timeBase = 1:45; winERP = 46:90;
        srcERP = zeros(numROIs, 90 );
        srcERP(ac_sources,winERP) = 1;        
        
        %% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj 
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));   
        Y_avg = Y;
        for iSub=1:numSubs
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
            sourceNoise = sourceData; % no signal, used to compute SNR
%             if length([idxROIfwd{iSub,ac_sources}]) ~= length(unique([idxROIfwd{iSub,ac_sources}]))
%                 fprintf('Overlapping source indexes S%d \n',listSub(iSub));
%             end
            for ss=1:length(ac_sources)
                % note that if there is overlapping index (same idx for 2
                % ROIs), the value in sourceData will be of the latest
                % source
                sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
            end
            % multiply fwd (128*20484) with the activated idx over time
            % (sourceData of 20484*90) and obtain Y elec x time
            y_stim = fullFwd{iSub} * sourceData;
%             % add noise
%             [noisy_data] = add_ERPnoise_with_SNR( y_stim , SNRlevel,winERP );
%             Y(iSub,:,:) = y_stim + noisy_data;
            Y(iSub,:,:) = y_stim;
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
        [betaLCFUN, betaCurv, betaBest, lambdaLCFUN, lambdaCurv, lambdaBest, ...
            lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(Y_avg,1)),srcERP);
        
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
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole,1));
        retrieveROI = squeeze(mean(regionROI,1));
        retrieveROIin = squeeze(mean(betaROIin,1));

        %% amount of crosstalk
        % leakage: amplitude signal in all areas, the "true" one is used
        % for normalisation (=1). Use RMS
        normTerm = rms(betaCurv(seedRoi,:));
        normTermBest = rms(betaBest(seedRoi,:));
        normTermWhole = rms(retrieveWhole(seedRoi,:));
        normTermROI = rms(retrieveROI(seedRoi,:));
        normTermROIin = rms(retrieveROIin(seedRoi,:));
        for iRoi = 1:numROIs
            crossTalkTemplate(repBoot,seedRoi,iRoi) = rms(betaCurv(iRoi,:)) / normTerm;
            crossTalkWhole(repBoot,seedRoi,iRoi) = rms(retrieveWhole(iRoi,:)) / normTermWhole;
            crossTalkROI(repBoot,seedRoi,iRoi) = rms(retrieveROI(iRoi,:)) / normTermROI;
            crossTalkROIin(repBoot,seedRoi,iRoi) = rms(retrieveROIin(iRoi,:)) / normTermROIin;
            crossTalkTemplateBest(repBoot,seedRoi,iRoi) = rms(betaBest(iRoi,:)) / normTermBest;
        end
             
        
    end % end boostrap
    
end % different activated sources

save('simulOutput/simulCrossTalk200.mat','crossTalkTemplate','crossTalkWhole',...
    'crossTalkROI','crossTalkROIin','crossTalkTemplateBest');

figure;
subplot(2,3,1);
imagesc(squeeze(mean(crossTalkTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))  
ylabel('seedArea');xlabel('predictArea')
title('templateBased')
subplot(2,3,2);imagesc(squeeze(mean(crossTalkWhole(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('Whole gcv')
subplot(2,3,3);imagesc(squeeze(mean(crossTalkROI(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('ROI gcv')
subplot(2,3,4);imagesc(squeeze(mean(crossTalkROIin(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('Oracle')
subplot(2,3,5);imagesc(squeeze(mean(crossTalkTemplateBest(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('template best reg param')
set(gcf,'position',[100,100,1800,800])
saveas(gcf,['figures' filesep 'crossTalkStepSNR' num2str(SNRlevel)],'png')

