clearvars;close all;
% same prog as simulV1MT but noise is added on forward models = correlated
% brain noise instead of adding it on the electrodes
% simulate ERP from V1, MT, V1+MT in 3 different time windows
% retrieve the sources using template method 
% Each simulation uses a different set of sbj and ERP but the same one is
% tested across different levels of noise and number of sbj
% Noise is simulated either on the electrodes (as earlier) or on the
% forward of each subject (which ends up with correlated noise = brain
% noise)

% 1-45 = baseline
% 46-90 = V1
% 91-135 = MT
% 136-180 = MT+V1

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

nbSbjToInclude =[2 8 20 50];

totBoot = 30;

for repBoot=1:totBoot
    
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
            % use the generated sources to simulate scalp activity for each sbj
            % (using individual fwd model)
            Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
            YbrainNoise = Y;
            Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
            Yb_avg = Y_avg;
            
            for iSub=1:numSubs
                % initialise matrix of source activity
                sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
                for ss=1:length(ac_sources)
                    % note that if there is overlapping index (same idx for 2
                    % ROIs), the value in sourceData will be of the latest
                    % source
                    sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
                end
                % multiply fwd (128*20484) with the activated idx over time
                % (sourceData of 20484*90) and obtain Y elec x time
                y_stim = fullFwd{iSub} * sourceData;
                % add noise to the fwd (= correlated brain noise)
                noiseGauss = 1/sqrt(2*noiseLevel) * randn(size(fullFwd{iSub},2),size(y_stim,2));
                brainNoise = fullFwd{iSub} * noiseGauss; 
                YbrainNoise(iSub,:,:) = y_stim + brainNoise;

                 % add noise to the electrodes
                [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
                % to keep the same SNR for the 2 Y, need to compute noise for
                % the 2 Y separately as it is based on the variance of the signal
                Y(iSub,:,:) = y_stim + noisy_data;
                
            end
           
            %%% Use average reference for centering Y
            %%% that is: substract the average electrode activity at each time point
            % this is done by bsxfun which applies element-wise substraction (the 90
            % averages across electrodes) - useless?
            for iSub=1:numSubs
                Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
                Yb_avg(iSub,:,:) = bsxfun(@minus,squeeze(YbrainNoise(iSub,:,:)), mean(squeeze(YbrainNoise(iSub,:,:))));
            end
 
            %% compute minimum norm
            regionWhole = zeros(numSubs,numROIs,length(srcERP));
            regionROI = regionWhole;
%             regionROILC = regionWhole;
%             regionWholeLC = regionWhole;
            betaROIin = regionWhole;
%             betaROIinLC = regionWhole;
            
            % min_norm on average data: get beta values for each ROI over time
            [betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
%             [betaAverageB, lambdaB] = minNormFast_lcurve(avMap, squeeze(mean(Yb_avg,1)));
            [betaLCFUN, betaAverageB, betaBest, lambdaB, lambdaCurv, lambdaBest, ...
                lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(Yb_avg,1)),srcERP);

%             indFwdROI_noise=[roiFwd{iSub,:}];
%             indData_noise=squeeze(Y_avg(iSub,:,:));
%             indFwd_noise=fullFwd{iSub};
            
            for iSub=1:numSubs
                % regular minimum_norm: on the 20484 indexes per sbj
                [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Yb_avg(iSub,:,:)), nLambdaRidge);
%                 [betaWholeLC, lambdaWholeLC] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
                
                [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Yb_avg(iSub,:,:)), nLambdaRidge);
%                 [betaROILC, lambdaROILC] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
                
                % beta values are for the indexes, but needs it per ROI
                % get the number of indexes per ROI for this subj
                rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
                % get the range
                range = [0 cumsum(rangeROI)]; % cumulative sum of elements
                % SUM (not average) the beta values per ROI (=across the indexes)
                regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%                 regionROILC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
                
                % need to find the indexes for whole brain -> use idxROIfwd
                % (no need to get the range)
                regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%                 regionWholeLC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
                
                % feed ROI per sbj instead of mesh = "oracle" (best possible recovery)
                sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
                [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Yb_avg(iSub,:,:)), nLambdaRidge);
%                 [betaROIinLC(iSub,:,:), lambdaGridMinNormROIinLC] = minNormFast_lcurve(sbjROI, squeeze(Y_avg(iSub,:,:)));
            end
            % average across subj
            retrieveWhole = squeeze(mean(regionWhole,1));
            retrieveROI = squeeze(mean(regionROI,1));
            retrieveROIin = squeeze(mean(betaROIin,1));
%             retrieveWholeLC = squeeze(mean(regionWholeLC,1));
%             retrieveROILC = squeeze(mean(regionROILC,1));
%             retrieveROIinLC = squeeze(mean(betaROIinLC,1));
            
            % save simulation
            simulERP(totSbj,level).listROIs = listROIs;
            simulERP(totSbj,level).listSub = listSub;
            simulERP(totSbj,level).winERP = winERP;
            simulERP(totSbj,level).srcERP = srcERP;
            simulERP(totSbj,level).data = Y_avg;
            simulERP(totSbj,level).dataCorrelNoise = Yb_avg;
            simulERP(totSbj,level).noise = SNRlevel(level);
            simulERP(totSbj,level).beta(1,:,:) = betaAverage;
            simulERP(totSbj,level).beta(2,:,:) = betaAverageB;
            simulERP(totSbj,level).beta(3,:,:) = retrieveWhole;
            simulERP(totSbj,level).beta(4,:,:) = retrieveROI;
            simulERP(totSbj,level).beta(5,:,:) = retrieveROIin;
            simulERP(totSbj,level).beta(6,:,:) = betaBest;
%             simulERP(repBoot,totSbj,level).beta(6,:,:) = retrieveWholeLC;
%             simulERP(repBoot,totSbj,level).beta(7,:,:) = retrieveROILC;
%             simulERP(repBoot,totSbj,level).beta(8,:,:) = retrieveROIinLC;
            
        end % noise
        
    end % sbj
    save(['simulOutput/brainNoise/simulV1MTbrainNoise' num2str(repBoot) '.mat'],simulERP)

end % boot




% nbModel = 2;
% % initialise variables
% aucAve = zeros(size(simulERP,1),size(simulERP,2),size(simulERP,3),nbModel);
% energyAve = aucAve;
% mseAveNorm = aucAve;
% 
% for model=1:nbModel
% for repBoot=1:size(simulERP,1)
%     for totSbj=1:size(simulERP,2)
%         for level=1:size(simulERP,3)            
%         [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
%             computeMetrics(squeeze(simulERP(repBoot,totSbj,level).beta(model,:,winERP)),simulERP(repBoot,totSbj,level).srcERP(:,winERP));        
%         end
%     end
% end
% end
% 
% %%% plot metrics
% modName = {'electrodeNoise','brainNoise'};
% figure;hold on
% for model=1:nbModel
% for ss=1:size(simulERP,2)
%     subplot(3,nbModel,model);hold on;
%     errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
%     xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
%     title(modName(model))
%     subplot(3,nbModel,model+nbModel);hold on;
%     errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
%     ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
%     subplot(3,nbModel,model+nbModel*2);hold on;
%     errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
%     ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
% end
% end
% legend('N=2','N=8','N=20','N=50')
% set(gcf,'position',[100 100 500 700])
% saveas(gcf,['figures' filesep 'brainNoiseV1MT'],'png')
