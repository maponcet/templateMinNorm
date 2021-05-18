clearvars;close all;
% simulation: V1+MT
% for a given simulation, amplitude and time function is different for each
% source but the same across participants (with different fwd models)
% simulation consistent with retinotopy: L&R sources are the same
% (= assumes full field stimulation)
% simulation on MESH
% minimum_norm: done on a)whole brain, b)only 18 ROIs, c)average50ROIs

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = 0.1; % noise level 10%
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm
numCols = 2; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

nbSbjToInclude =[1 2 4 8 16 25 32 50];
totBoot = 20; % nb of bootstrap

% initialise variables
aucMinNorm = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
aucGpMinNorm = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
aucMinNormInd = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
energyMinNorm = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
energyGpNorm = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
energyMinNormInd = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
mseMinNorm = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
mseGpNorm = zeros(totBoot,length(nbSbjToInclude),length(activeROIs));
mseMinNormInd= zeros(totBoot,length(nbSbjToInclude),length(activeROIs));


for totSbj=1:length(nbSbjToInclude)
    numSubs = nbSbjToInclude(totSbj);
    
    for repBoot =1:totBoot
        fprintf('sbj%d bootstrap %d \n',numSubs,repBoot);
        clear Y Ylo Yerp YloERP
        
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
        % amplitude and time function is different for each source but
        % the same for all sbj for a given bootstrap
        [sourceAmplitude, sourceERP] = createSourceROI(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
        
        %% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj 
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),length(sourceERP));
        for iSub=1:numSubs
            % initialise matrix of source activity
            sourceData = zeros(size(fullFwd{iSub},2) , length(sourceERP));
            if length([idxROIfwd{iSub,ac_sources}]) ~= length(unique([idxROIfwd{iSub,ac_sources}]))
                fprintf('Overlapping source indexes S%d \n',listSub(iSub));
            end
            for ss=1:length(ac_sources)
                % note that if there is overlapping index (same idx for 2
                % ROIs), the value in sourceData will be of the latest
                % source
                sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(sourceERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
            end
            % multiply fwd (128*20484) with the activated idx over time
            % (sourceData of 20484*90) and obtain Y elec x time
            y_stim = fullFwd{iSub} * sourceData;
            % add noise
            [noisy_data,~] = add_noise_with_SNR( y_stim , SNRlevel );
            Y(iSub,:,:) = y_stim + noisy_data;
        end
        %         figure;plotContourOnScalp(squeeze(Y(1,:,45)),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
        
        %% center Y, fullFwd, roiFwd across electrodes?? Sbj??
        %%%%%% what should be substracted??????
        % reshape Y?
        % for now centre for each sbj across electrodes
        % IT IS USELESS TO DO IT FOR fullFwd or roiFwd if it is only on
        % single sbj -> gives the same numbers!!!
        for iSub=1:numSubs
            Y(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
        end
        %             fullFwd{iSub} = bsxfun(@minus,fullFwd{iSub}, mean(fullFwd{iSub}));
        %             fullFwd2 = bsxfun(@minus,[fullFwd{:}], mean([fullFwd{:}]));
        %             for rr=1:numROIs
        %                 roiFwd{iSub,rr} = bsxfun(@minus,roiFwd{iSub,rr}, mean(roiFwd{iSub,rr}));
        %             end
        
        %%  reduce dimension data (Ylo)
        % =PCA denoised version of Y (denoised by truncation of the SVD)
        %%%%%% per sbj or alltogether?
        for iSub=1:numSubs
            [u1, s1, v1] = svd(squeeze(Y(iSub,:,:)));
            Ylo(iSub,:,:) = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
        end
        %         figure;plotContourOnScalp(YloTEST(:,45),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
        %         figure;plotContourOnScalp(y_stim(:,45),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
       
        
        
        %% compute minimum norm
        regionWhole = zeros(numSubs,numROIs,length(sourceERP));
        regionROI = zeros(numSubs,numROIs,length(sourceERP));        
        % min_norm on average data: get beta values for each ROI over time 
        [betaAverage, ~, lambdaAverage] = minimum_norm(avMap, squeeze(mean(Ylo)), nLambdaRidge);
        for iSub=1:numSubs
            % regular minimum_norm: on the 20484 indexes per sbj
            [betaWhole, ~, lambdaWhole] = minimum_norm(fullFwd{iSub}, squeeze(Ylo(iSub,:,:)), nLambdaRidge);
            % min_norm on only the ROI indexes per sbj
            [betaROI, ~, lambdaROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Ylo(iSub,:,:)), nLambdaRidge);
            
            % beta values are for the indexes, but I want it per ROI
            % get the number of indexes per ROI for this subj
            rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
            % get the range
            range = [0 cumsum(rangeROI)]; % cumulative sum of elements
            % average the beta values per ROI (=across the indexes)
            regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); 
            % need to find the indexes for whole brain
            regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
            % create brain resp for each sbj
            % YhatWhole(iSub,:,:) = fullFwd{iSub}*betaWhole(iSub,:,:);
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole));
        retrieveROI = squeeze(mean(regionROI));
        
        %% plot source amplitude overtime per ROI
        % for 1)source signal 2)trad+ROI informed average whole brain 
        % 3) trad ROI 4)new average-based
        %%%%%%%%%%%%% plot simulated signal instead of source!!!! 
        %%%%%%% but for this I would need to get it when simulating Y and
        %%%%%%% average the ROI indexes
        if numSubs == 50 && repBoot == 1
            figure;rr=1;
            set(gcf,'position',[100,100,800,1000])
            for iRoi = 1:2:numROIs
                subplot(9,4,1+4*(rr-1));hold on;
                plot(sourceERP(iRoi,:)','m-','linewidth',2)
                plot(sourceERP(iRoi+1,:)','b-','linewidth',2)
                ylim([min(min(sourceERP)) max(max(sourceERP))])
                title(listROIs(iRoi))
                subplot(9,4,2+4*(rr-1));hold on;
                plot(retrieveWhole(iRoi,:)','m-','linewidth',2)
                plot(retrieveWhole(iRoi+1,:)','b-','linewidth',2)
                ylim([min(min(retrieveWhole)) max(max(retrieveWhole))])
                subplot(9,4,3+4*(rr-1));hold on;
                plot(retrieveROI(iRoi,:)','m-','linewidth',2)
                plot(retrieveROI(iRoi+1,:)','b-','linewidth',2)
                ylim([min(min(retrieveROI)) max(max(retrieveROI))])
                subplot(9,4,4+4*(rr-1));hold on;
                plot(betaAverage(iRoi,:)','m-','linewidth',2)
                plot(betaAverage(iRoi+1,:)','b-','linewidth',2)
                ylim([min(min(betaAverage)) max(max(betaAverage))])                
                rr=rr+1;
            end
            saveas(gcf,['figures' filesep 'simulMay'],'png')
        end       
        
        
        %% compute auc, mse, relative energy using average signal in rois
        % do for all the min norm outputs
        [aucAve(repBoot,totSbj), energyAve(repBoot,totSbj),...
            mseAve(repBoot,totSbj)] = computeMetrics(betaAverage,ac_sources,sourceERP);
        [aucWhole(repBoot,totSbj), energyWhole(repBoot,totSbj),...
            mseWhole(repBoot,totSbj)] = computeMetrics(retrieveWhole,ac_sources,sourceERP);
        [aucROI(repBoot,totSbj), energyROI(repBoot,totSbj),...
            mseROI(repBoot,totSbj)] = computeMetrics(retrieveROI,ac_sources,sourceERP);

        
    end % end boostrap
    
end % go through nb of sub to include



nameSim = {'V1-MT','V2V-V4'};
for simul=1:length(allSig)
    figure;
    subplot(1,3,1+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(aucOld(:,:,simul)),std(aucOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucNew(:,:,simul)),std(aucNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucMinNormInd(:,:,simul)),std(aucMinNormInd(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucOldERP(:,:,simul)),std(aucOldERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucNewERP(:,:,simul)),std(aucNewERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('AUC')
    ylim([0 1])
    title(['simul ' nameSim{simul}])
    subplot(1,3,2+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(energyOld(:,:,simul)),std(energyOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyNew(:,:,simul)),std(energyNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyMinNormInd(:,:,simul)),std(energyMinNormInd(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyOldERP(:,:,simul)),std(energyOldERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyNewERP(:,:,simul)),std(energyNewERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('Energy')
    ylim([0 1])
    subplot(1,3,3+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(mseOld(:,:,simul)),std(mseOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseNew(:,:,simul)),std(mseNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseMinNormInd(:,:,simul)),std(mseMinNormInd(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseOldERP(:,:,simul)),std(mseOldERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseNewERP(:,:,simul)),std(mseNewERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('MSE')
    ylim([0 2])
    legend('minNorm','newMethod','minNormInd','minNormERP','newMethodERP')
end
set(gcf,'position',[100,100,900,500])
saveas(gcf,['figures' filesep 'testSimulMetrics_Retino-new'],'png')

