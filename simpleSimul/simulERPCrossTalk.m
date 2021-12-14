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
SNRlevel = 10; % noise level, 10 times more signal than noise
totBoot = 10; % nb of bootstrap
crossTalkTemplate = zeros(totBoot,numROIs,numROIs);
crossTalkWhole = zeros(totBoot,numROIs,numROIs);
crossTalkROI = zeros(totBoot,numROIs,numROIs);

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
        
        
        %% Simulate sources (sourceERP)
        % amplitude (1 to 10) and time function is different for each 
        % source but the same for all sbj for a given bootstrap
        [srcAmp, srcSSVEP, srcERP,winERP] = createSourceROI(numROIs,ac_sources,[]);
        % ERP baseline timewindow
        timeBase = setdiff(1:size(srcERP,2),winERP);
        
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
        regionWhole = zeros(numSubs,numROIs,length(srcERP));
        regionROI = zeros(numSubs,numROIs,length(srcERP));        
        % min_norm: get beta values for each ROI over time 
        % stack all the participants vertically
        stackY = reshape(permute(Y_avg,[2,1,3]),[size(Y_avg,1)*size(Y_avg,2),size(Y_avg,3)]);
        % stack the template for as many participants
        stackAvMap = repmat(avMap,numSubs,1);
        [betaAverage, betaMinNorm, lambda, gcvErrorMinNorm, lambdaGridMinNorm] = minimum_norm(stackAvMap, stackY, nLambdaRidge);

        for iSub=1:numSubs
            clear betaWhole betaROI
            % regular minimum_norm: on the 20484 indexes per sbj
            [betaWhole, ~, lambdaWhole] = minimum_norm(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            % min_norm on only the ROI indexes per sbj
            [betaROI, ~, lambdaROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);

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
                        
%             % check ind topo....
%             figure;subplot(1,3,1)
%             plotOnEgi(squeeze(Y(iSub,:,nT)))
%             subplot(1,3,2)
%             yWhole(iSub,:,:) = [fullFwd{iSub}] * betaWhole;
%             plotOnEgi(squeeze(yWhole(iSub,:,nT)))
%             subplot(1,3,3)
%             yROI(iSub,:,:) = [roiFwd{iSub,:}] * betaROI;
%             plotOnEgi(squeeze(yROI(iSub,:,nT)))
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole,1));
        retrieveROI = squeeze(mean(regionROI,1));

        
        %% amount of crosstalk
        % leakage: amplitude signal in all areas, the "true" one is used
        % for normalisation (=1). Use RMS
        normTerm = rms(betaAverage(seedRoi,:));
        normTermWhole = rms(retrieveWhole(seedRoi,:));
        normTermROI = rms(retrieveROI(seedRoi,:));
        for iRoi = 1:numROIs
%             crossTalkTemplate(repBoot,seedRoi,iRoi) = rms(betaAverage(iRoi,:)) / normTerm;
%             crossTalkWhole(repBoot,seedRoi,iRoi) = rms(retrieveWhole(iRoi,:)) / normTermWhole;
%             crossTalkROI(repBoot,seedRoi,iRoi) = rms(retrieveROI(iRoi,:)) / normTermROI;
            test(seedRoi,iRoi) = rms(betaAverage(iRoi,:)) / normTerm;
        end
            
%         count = 1;
%         figure;set(gcf,'position',[100,100,800,1000])
%         for iRoi = 1:2:numROIs
%             % need to normalise the signal
%             subplot(3,3,count);hold on
%             plot(srcERP(iRoi,:) / max(max(abs(srcERP))) ,'LineWidth',2);
%             plot(srcERP(iRoi+1,:) / max(max(abs(srcERP))) ,'LineWidth',2);
%             tt = cell2mat(listROIs(iRoi));title(tt(1:end-2) ,'LineWidth',2)
%             ylim([-1 1]);count=count+1;
%         end
%         count = 1;
%         figure;set(gcf,'position',[100,100,800,1000])
%         for iRoi = 1:2:numROIs
%             % need to normalise the signal
%             subplot(3,3,count);hold on
%             plot(betaAverage(iRoi,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
%             plot(betaAverage(iRoi+1,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
%             tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%             ylim([-1 1]);count=count+1;
%         end        
%         count = 1;
%         figure;set(gcf,'position',[100,100,800,1000])
%         for iRoi = 1:2:numROIs
%             % need to normalise the signal
%             subplot(3,3,count);hold on
%             plot(retrieveROI(iRoi,:) / max(max(abs(retrieveROI))) ,'LineWidth',2);
%             plot(retrieveROI(iRoi+1,:) / max(max(abs(retrieveROI))) ,'LineWidth',2);
%             tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%             ylim([-1 1]);count=count+1;
%         end   
        
    end % end boostrap
    
end % different activated sources

figure;
subplot(3,1,1);
imagesc(squeeze(mean(crossTalkTemplate(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))  
ylabel('seedArea');xlabel('predictArea')
title('templateBased')
subplot(3,1,2);imagesc(squeeze(mean(crossTalkWhole(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('Whole MinNorm')
subplot(3,1,3);imagesc(squeeze(mean(crossTalkROI(:,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
ylabel('seedArea');xlabel('predictArea')
title('ROI MinNorm')
set(gcf,'position',[100,100,800,2000])
saveas(gcf,['figures' filesep 'ERPcrossTalkSNR' num2str(SNRlevel)],'png')


% figure;
% subplot(3,1,1);imagesc(squeeze(crossTalkTemplate(1,[1:2:18 2:2:18],[1:2:18 2:2:18])));colorbar;%caxis([0 1])
% set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
% set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))  
% ylabel('seedArea');xlabel('predictArea')
% title('templateBased')
% subplot(3,1,2);imagesc(squeeze((crossTalkWhole(1,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
% set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
% set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
% ylabel('seedArea');xlabel('predictArea')
% title('Whole MinNorm')
% subplot(3,1,3);imagesc(squeeze((crossTalkROI(1,[1:2:18 2:2:18],[1:2:18 2:2:18]))));colorbar;%caxis([0 1])
% set(gca, 'XTick',1:18, 'XTickLabel',listROIs([1:2:18 2:2:18]))   
% set(gca, 'YTick',1:18, 'YTickLabel',listROIs([1:2:18 2:2:18]))   
% ylabel('seedArea');xlabel('predictArea')
% title('ROI MinNorm')
% set(gcf,'position',[100,100,800,1800])
% saveas(gcf,['figures' filesep 'ERPcrossTalkNoNoise2'],'png')