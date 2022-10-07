clearvars;close all;
% Test for an interfering source not in visual areas - picked randomly from
% anatomical ROIs
% ERP activity from V1+MT 
% scaling of the external ROI by 50%: take the left external ROI / MT-L
% always N=50

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128_allROI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
load('nonVisualROI.mat') % load list of non visual ROIs

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};

% sourceExt = {'caudalanteriorcingulate-L'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

numSubs = 50;
totBoot = 60;

for repBoot=1:totBoot
    % select external source
    srcExtInd = randi(length(nonVisualROI)/2,1,1);
    sourceExt = [nonVisualROI(srcExtInd) nonVisualROI(srcExtInd+length(nonVisualROI)/2)];
    
    % Simulate sources
    % ERP & baseline timewindow
    timeBase = 1:45;
    winERP = 46:90;
    srcERP = zeros(numROIs,45*2); % 45*2 timepoints
    srcExtERP = zeros(length(sourceExt),45*2);
    % create 1 ERP, repeated for the 2nd source (bilateral)
    srcERP(:,46:90) = createSourceERP(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
    srcExtERP(:,46:90) = createSourceERP(length(sourceExt),1,2); 

            
%     for totSbj=1:length(nbSbjToInclude) 
%         numSubs = nbSbjToInclude(totSbj);
        fprintf('N%d bootstrap %d\n',numSubs,repBoot)
        
        % list of random sbj with replacement
        listSub = randi(length(dirList),numSubs,1);
        
        %% LOAD FWD
        fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);normFwd=cell(1,numSubs);normFwd50=cell(1,numSubs);
        idxROIfwd=cell(numSubs,numROIs);idxExtROIfwd = cell(numSubs,length(sourceExt));
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
            % save index for interfering/external source
            for extSrc=1:length(sourceExt)
                indexExtROI = find(strcmp(sourceExt(extSrc),{roiInfo.name}));
                idxExtROIfwd{iSub,extSrc} = roiInfo(indexExtROI).meshIndices;
            end

            % normalise forward for equal source strengh
            roiNormVis = vecnorm(fwdMatrix(:,idxROIfwd{iSub,17})'); % MT-L
            roiNormExt = vecnorm(fwdMatrix(:,idxExtROIfwd{iSub,1})'); % first interfering ROI
%             figure;plot(roiNormVis);hold on; plot(roiNorExt)
            scalingFactor = roiNormVis / roiNormExt; % mrdivide
%             fwdMatrixN = fwdMatrix; 
%             fwdMatrixN(:,idxExtROIfwd{iSub,1}) = scalingFactor * fwdMatrix(:,idxExtROIfwd{iSub,1}); 
%             normFwd{iSub} = fwdMatrixN;
            fwdMatrixN2 = fwdMatrix; 
            fwdMatrixN2(:,idxExtROIfwd{iSub,1}) = scalingFactor/2 * fwdMatrix(:,idxExtROIfwd{iSub,1}); 
            fwdMatrixN2(:,idxExtROIfwd{iSub,2}) = scalingFactor/2 * fwdMatrix(:,idxExtROIfwd{iSub,2}); 
            normFwd50{iSub} = fwdMatrixN2;
        end

        
        for level=1:length(SNRlevel)
            noiseLevel = SNRlevel(level);
            
            %%% Simulate scalp activity (Y)
            % use the generated sources to simulate scalp activity for each sbj
            % (using individual fwd model)
            Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
            Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
            Yext = Y; Y_avgExt = Y_avg;Yext50 = Y; Y_avgExt50 = Y_avg;
            
            for iSub=1:numSubs
                % initialise matrix of source activity
                sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
                
                for ss=1:length(ac_sources)
                    % note that if there is overlapping index (same idx for 2
                    % ROIs), the value in sourceData will be of the latest
                    % source
                    sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
                end
                % add interfering source to the sourceData
                 sourceDataExt = sourceData; 
               for ss=1:length(sourceExt)
                    sourceDataExt(idxExtROIfwd{iSub,ss},:) = repmat(srcExtERP(ss,:),length(idxExtROIfwd{iSub,ss}),1);
                end
                % multiply fwd (128*20484) with the activated idx over time
                % (sourceData of 20484*90) and obtain Y elec x time
                y_stim = fullFwd{iSub} * sourceData;
%                 y_stimExt = normFwd{iSub} * sourceDataExt;
                y_stimExt50 = normFwd50{iSub} * sourceDataExt;
                % add gaussian noise on electrodes
                [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
%                 [noisyData] = add_ERPnoise_with_SNR( y_stimExt , noiseLevel,winERP );
                [noisyData2] = add_ERPnoise_with_SNR( y_stimExt50 , noiseLevel,winERP );
                % to keep the same SNR for the 2 Y, need to compute noise for
                % the 2 Y separately as it is based on the variance of the signal
                Y(iSub,:,:) = y_stim + noisy_data;
%                 Yext(iSub,:,:) = y_stimExt + noisyData;
                Yext50(iSub,:,:) = y_stimExt50 + noisyData2;

%             % check SNR
%             allSig = Y(iSub,:,:); erpSig = Y(iSub,:,winERP); erpNoise = noisy_data(:,winERP);
%             (rms(allSig(:))/rms(noisy_data(:)))^2 - 1 % half: 1/2win=no sig
%             (rms(erpSig(:))/rms(erpNoise(:)))^2 - 1 % should = SNR
%             allSig = Yext50(iSub,:,:); erpSig = Yext50(iSub,:,winERP); erpNoise = noisyData2(:,winERP);
%             (rms(allSig(:))/rms(noisyData2(:)))^2 - 1 % half: 1/2win=no sig
%             (rms(erpSig(:))/rms(erpNoise(:)))^2 - 1 % should = SNR
            end
            
            %%% Use average reference for centering Y
            %%% that is: substract the average electrode activity at each time point
            % this is done by bsxfun which applies element-wise substraction (the 90
            % averages across electrodes) - useless?
            for iSub=1:numSubs
                Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
%                 Y_avgExt(iSub,:,:) = bsxfun(@minus,squeeze(Yext(iSub,:,:)), mean(squeeze(Yext(iSub,:,:))));
                Y_avgExt50(iSub,:,:) = bsxfun(@minus,squeeze(Yext50(iSub,:,:)), mean(squeeze(Yext50(iSub,:,:))));
            end
 
            %% compute minimum norm
            regionWhole = zeros(numSubs,numROIs,length(srcERP));
            regionROI = regionWhole;
            betaROIin = regionWhole;
%             regionROILC = regionWhole;
%             regionWholeLC = regionWhole;
%             betaROIinLC = regionWhole;
            
            % min_norm on average data: get beta values for each ROI over time
            [betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
%             [betaAverageExt, lambdaExt] = minNormFast_lcurve(avMap, squeeze(mean(Y_avgExt,1)));
            [betaAverageExt50, lambdaExt50] = minNormFast_lcurve(avMap, squeeze(mean(Y_avgExt50,1)));
            
% %             indFwdROI_noise=[roiFwd{iSub,:}];
% %             indData_noise=squeeze(Y_avg(iSub,:,:));
% %             indFwd_noise=fullFwd{iSub};
%             
%             for iSub=1:numSubs
%                 % regular minimum_norm: on the 20484 indexes per sbj
%                 [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avgExt50(iSub,:,:)), nLambdaRidge);
% %                 [betaWholeLC, lambdaWholeLC] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
%                 
%                 [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avgExt50(iSub,:,:)), nLambdaRidge);
% %                 [betaROILC, lambdaROILC] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
%                 
%                 % beta values are for the indexes, but needs it per ROI
%                 % get the number of indexes per ROI for this subj
%                 rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
%                 % get the range
%                 range = [0 cumsum(rangeROI)]; % cumulative sum of elements
%                 % SUM (not average) the beta values per ROI (=across the indexes)
%                 regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
% %                 regionROILC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%                 
%                 % need to find the indexes for whole brain -> use idxROIfwd
%                 % (no need to get the range)
%                 regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
% %                 regionWholeLC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%                 
%                 % feed ROI per sbj instead of mesh = "oracle" (best possible recovery)
%                 sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
%                 [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avgExt50(iSub,:,:)), nLambdaRidge);
% %                 [betaROIinLC(iSub,:,:), lambdaGridMinNormROIinLC] = minNormFast_lcurve(sbjROI, squeeze(Y_avg(iSub,:,:)));
%             end
%             % average across subj
%             retrieveWhole = squeeze(mean(regionWhole,1));
%             retrieveROI = squeeze(mean(regionROI,1));
%             retrieveROIin = squeeze(mean(betaROIin,1));
% %             retrieveWholeLC = squeeze(mean(regionWholeLC,1));
% %             retrieveROILC = squeeze(mean(regionROILC,1));
% %             retrieveROIinLC = squeeze(mean(betaROIinLC,1));
            
            % save simulation
            simulERP(level).nonVisualROI = sourceExt(1);
            simulERP(level).listROIs = listROIs;
            simulERP(level).listSub = listSub;
            simulERP(level).winERP = winERP;
            simulERP(level).srcERP = srcERP;
            simulERP(level).data = Y_avg;
            simulERP(level).dataExt = Y_avgExt50;
            simulERP(level).noise = SNRlevel(level);
            simulERP(level).beta(1,:,:) = betaAverage;
%             simulERP(repBoot,totSbj,level).beta(2,:,:) = betaAverageExt;
            simulERP(level).beta(2,:,:) = betaAverageExt50;
%             simulERP(repBoot,totSbj,level).beta(4,:,:) = retrieveWhole;
%             simulERP(repBoot,totSbj,level).beta(5,:,:) = retrieveROI;
%             simulERP(repBoot,totSbj,level).beta(6,:,:) = retrieveROIin;
%             simulERP(repBoot,totSbj,level).beta(5,:,:) = retrieveWholeLC;
%             simulERP(repBoot,totSbj,level).beta(6,:,:) = retrieveROILC;
%             simulERP(repBoot,totSbj,level).beta(7,:,:) = retrieveROIinLC;
            
        end % noise
        
%     end % sbj

    save(['simulOutput/interference/simulV1MTinterferenceERP' num2str(repBoot) '.mat'],'simulERP')
    
end % boot

