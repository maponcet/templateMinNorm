clearvars;close all;
% simulation1: V1+MT
% simulation2: V2v+V4
% fix amplitude and waveform for a given simulation (difference is the
% participant's fwd) but they randomly change across bootstrap
% UPDATE simulation should be consistent with retinotopy: L&R sources 
% should be the same (assuming full field stimulation) as well as ventral
% and dorsal
% This program only includes Left areas (assumes right stim). Retrieval is
% done on either all areas or only left areas

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50.mat') % load average map of ROIs (128 elec x 18 ROIs)
listROIs = avMap.roiNames; % the order of the list MUST be the same as in the avMap!
numROIs = length(listROIs);

% some parameters
SNRlevel = 0.1; % noise level 10%
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm
numCols = 2; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
allSig = {{'V1-L','MT-L'}}; 
nbSbjToInclude =[1 2 4 8 16 25 32 50];
totBoot = 20; % nb of bootstrap

% initialise variables
aucNew = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
aucOld = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyOld = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyNew = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseOld = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseNew = zeros(totBoot,length(nbSbjToInclude),length(allSig));

aucNewUni = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
aucOldUni = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyOldUni = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyNewUni = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseOldUni = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseNewUni = zeros(totBoot,length(nbSbjToInclude),length(allSig));   

aucNewERP = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
aucOldERP = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyOldERP = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyNewERP = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseOldERP = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseNewERP = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 

for simulation = 1:length(allSig)
    % simulated signal
    signalROIs = allSig{simulation};
    
    % find the ROI index which corresponds to the one in signalROIs
    ac_sources = cell2mat(arrayfun(@(x) cellfind(avMap.roiNames,signalROIs{x}),1:length(signalROIs),'uni',false));
    
    for totSbj=1:length(nbSbjToInclude)
        numSubs = nbSbjToInclude(totSbj);
        
        for repBoot =1:totBoot
            fprintf('sim%d sbj%d bootstrap %d \n',simulation,numSubs,repBoot);
            clear Y Ylo Yerp YloERP
            
            % list of random sbj with replacement
            listSub = randi(length(dirList),numSubs,1);
            
            %% GET ROIs
            % do for each sbj separately but on 18 ROIs, not mesh
            % so first need to get the ROI per sbj
            roiMap=zeros(128,numROIs,numSubs);
            roiStacked = [];
            for iSubj=1:numSubs
                clear fwdMatrix roiInfo
                % fwd file
                load([dataPath dirList(listSub(iSubj)).name])
                % check that all rois are there
                if sum(ismember(listROIs, {roiInfo.name})) ~= numROIs
                    fprintf('nb of ROIs mismatch for %s \n',listFiles(iSubj).name)
                else
                    % go through each ROI and average all the mesh indexes corresponding to that ROI
                    for rr=1:numROIs
                        clear indexROI
                        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
                        [roiMap(:,rr,iSubj)] = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
                    end
                    roiStacked = blkdiag(roiStacked, roiMap(:,:,iSubj));
                end
            end
            roiStacked = bsxfun(@minus,roiStacked, mean(roiStacked));
            
            %% Simulate sources
            % amplitude and time function is different for each source but
            % the same across participants for a given bootstrap
            [Y, Yerp, sourceValOverTime, sourceValOverTimeERP] = simulSource(roiMap,SNRlevel, ac_sources,[]);
            
            %%  reduce dimension data (Ylo)
            % =PCA denoised version of Y (denoised by truncation of the SVD)
            n = numel(Y);
            [u1, s1, v1] = svd(Y);
            Ylo = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
            % compute average EEG 
            unstackedData = reshape(Ylo,128,numSubs,size(Ylo,2));
            grandMeanData = squeeze(mean(unstackedData,2));
            
            [u1ERP, s1ERP, v1ERP] = svd(Yerp);
            YloERP = u1ERP(:,1:numCols)*s1ERP(1:numCols,1:numCols)*v1ERP(:, 1:numCols)';
            % compute average EEG 
            unstackedDataERP = reshape(YloERP,128,numSubs,size(YloERP,2));
            grandMeanDataERP = squeeze(mean(unstackedDataERP,2));
           
            %% compute minimum norm
            % retrieve the sources from all areas 
            [betaReg, ~, ~] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
            [betaRegERP, ~, ~] = minimum_norm(avMap.activity, grandMeanDataERP, nLambdaRidge);
            [betaMinNorm, ~, ~] = minimum_norm(roiStacked, Ylo, nLambdaRidge);
            [betaMinNormERP, ~, ~] = minimum_norm(roiStacked, YloERP, nLambdaRidge);
            % or only from one hemisphere
            % find the index for the left regions
            indexLeft = find(contains(avMap.roiNames,'-L'));
            [betaRegUni, ~, ~] = minimum_norm(avMap.activity(:,indexLeft), grandMeanData, nLambdaRidge);
            [betaRegERPUni, ~, ~] = minimum_norm(avMap.activity(:,indexLeft), grandMeanDataERP, nLambdaRidge);
            % indexes for the left regions for all stacked sbj
            allLeftIdx=[];
            for iSbj=1:numSubs
                allLeftIdx = [allLeftIdx indexLeft+18*(iSbj-1)];
            end
            [betaMinNormUni, ~, ~] = minimum_norm(roiStacked(:,allLeftIdx), Ylo, nLambdaRidge);
            [betaMinNormERPUni, ~, ~] = minimum_norm(roiStacked(:,allLeftIdx), YloERP, nLambdaRidge);
            
            % betaMinNorm is numSubs*numROIs x time so need to extract betaMinNorm of a single sbj
            % for each ROI separately then average across sbj
            regionTmp = arrayfun(@(x) betaMinNorm(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
            regionMinNorm =reshape(cell2mat(regionTmp),numROIs,size(Ylo,2),numSubs);
            regionMinNormAvg = mean(regionMinNorm,3);
            regionTmpERP = arrayfun(@(x) betaMinNormERP(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
            regionMinNormERP =reshape(cell2mat(regionTmpERP),numROIs,size(Ylo,2),numSubs);
            regionMinNormAvgERP = mean(regionMinNormERP,3);     
            % Divide numROIs by 2 since I am only working with half ROIs
            % for Uni
            regionTmpUni = arrayfun(@(x) betaMinNormUni(1+numROIs/2*(x-1):numROIs/2*x,:),1:numSubs,'uni',false);
            regionMinNormUni =reshape(cell2mat(regionTmpUni),numROIs/2,size(Ylo,2),numSubs);
            regionMinNormAvgUni = mean(regionMinNormUni,3);
            regionTmpERPUni = arrayfun(@(x) betaMinNormERPUni(1+numROIs/2*(x-1):numROIs/2*x,:),1:numSubs,'uni',false);
            regionMinNormERPUni =reshape(cell2mat(regionTmpERPUni),numROIs/2,size(Ylo,2),numSubs);
            regionMinNormAvgERPUni = mean(regionMinNormERPUni,3); 
            
            %% compute auc, mse, relative energy using average signal in rois
            % do for all the min norm outputs
            [aucNew(repBoot,totSbj,simulation), energyNew(repBoot,totSbj,simulation),...
                mseNew(repBoot,totSbj,simulation)] = computeMetrics(betaReg,ac_sources,sourceValOverTime);
            [aucNewERP(repBoot,totSbj,simulation), energyNewERP(repBoot,totSbj,simulation),...
                mseNewERP(repBoot,totSbj,simulation)] = computeMetrics(betaRegERP,ac_sources,sourceValOverTimeERP);
            [aucOld(repBoot,totSbj,simulation), energyOld(repBoot,totSbj,simulation),...
                mseOld(repBoot,totSbj,simulation)] = computeMetrics(regionMinNormAvg,ac_sources,sourceValOverTime);
            [aucOldERP(repBoot,totSbj,simulation), energyOldERP(repBoot,totSbj,simulation),...
                mseOldERP(repBoot,totSbj,simulation)] = computeMetrics(regionMinNormAvgERP,ac_sources,sourceValOverTimeERP);
            % index of sources when only 9 are considered
            idxSources = find(ismember(indexLeft,ac_sources));
            [aucNewUni(repBoot,totSbj,simulation), energyNewUni(repBoot,totSbj,simulation),...
                mseNewUni(repBoot,totSbj,simulation)] = computeMetrics(betaRegUni,idxSources,sourceValOverTime(indexLeft,:));        
            [aucOldUni(repBoot,totSbj,simulation), energyOldUni(repBoot,totSbj,simulation),...
                mseOldUni(repBoot,totSbj,simulation)] = computeMetrics(regionMinNormAvgUni,idxSources,sourceValOverTime(indexLeft,:));

            
            %% plot source amplitude overtime per ROI
            % for 1)true signal 2)new template minnorm 3)standard minnorm
            if numSubs == 50 && repBoot == 1
                figure;rr=1;
                set(gcf,'position',[100,100,1000,1000])
                for iRoi = 1:2:size(betaReg,1)
                    subplot(9,6,1+6*(rr-1));hold on;
                    plot(sourceValOverTime(iRoi,:)','m-','linewidth',2)
                    plot(sourceValOverTime(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(sourceValOverTime)) max(max(sourceValOverTime))])
                    title(listROIs(iRoi))                    
                    subplot(9,6,2+6*(rr-1));hold on;
                    plot(betaReg(iRoi,:)','m-','linewidth',2)
                    plot(betaReg(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(betaReg)) max(max(betaReg))])
                    subplot(9,6,3+6*(rr-1));hold on;
                    plot(regionMinNormAvg(iRoi,:)','m-','linewidth',2)
                    plot(regionMinNormAvg(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(regionMinNormAvg)) max(max(regionMinNormAvg))])
                    subplot(9,6,4+6*(rr-1));hold on;
                    plot(sourceValOverTimeERP(iRoi,:)','m-','linewidth',2)
                    plot(sourceValOverTimeERP(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(sourceValOverTimeERP)) max(max(sourceValOverTimeERP))])
                    subplot(9,6,5+6*(rr-1));hold on;
                    plot(betaRegERP(iRoi,:)','m-','linewidth',2)
                    plot(betaRegERP(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(betaRegERP)) max(max(betaRegERP))])
                    subplot(9,6,6+6*(rr-1));hold on;
                    plot(regionMinNormAvgERP(iRoi,:)','m-','linewidth',2)
                    plot(regionMinNormAvgERP(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(regionMinNormAvgERP)) max(max(regionMinNormAvgERP))])
                    rr=rr+1;
                end
                saveas(gcf,['figures' filesep 'testSimul_UniAll' num2str(simulation)],'png')
                figure;rr=1;
                set(gcf,'position',[100,100,800,1000])
                for iRoi = 1:size(betaReg,1)/2
                    subplot(3,3,rr);hold on;
                    plot(sourceValOverTimeERP(iRoi+1*(rr-1),:)','g-','linewidth',2)
                    ylim([min(min(sourceValOverTimeERP)) max(max(sourceValOverTimeERP))])
                    title(listROIs(iRoi+1*(rr-1)))
                    plot(regionMinNormAvgERPUni(iRoi,:)','b-','linewidth',2)
                    plot(betaRegERPUni(iRoi,:)','r-','linewidth',2)
                    rr=rr+1;
                end
                legend('source','minNorm','new','location','best')
                saveas(gcf,['figures' filesep 'testSimulERP_Uni' num2str(simulation)],'png')
            end
                       
        end % end boostrap
        
    end % go through numSub
    
end % different activated sources


nameSim = {'V1-MT','V2V-V4'};
for simul=1:length(allSig)
    figure;
    subplot(1,3,1+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(aucOld(:,:,simul)),std(aucOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucNew(:,:,simul)),std(aucNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucOldERP(:,:,simul)),std(aucOldERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucNewERP(:,:,simul)),std(aucNewERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('AUC')
    ylim([0 1])
    title(['simul ' nameSim{simul}])
    subplot(1,3,2+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(energyOld(:,:,simul)),std(energyOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyNew(:,:,simul)),std(energyNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyOldERP(:,:,simul)),std(energyOldERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyNewERP(:,:,simul)),std(energyNewERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('Energy')
    ylim([0 1])
    subplot(1,3,3+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(mseOld(:,:,simul)),std(mseOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseNew(:,:,simul)),std(mseNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseOldERP(:,:,simul)),std(mseOldERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseNewERP(:,:,simul)),std(mseNewERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('MSE')
    ylim([0 2])
    legend('minNorm','newMethod','minNormERP','newMethodERP')
end
set(gcf,'position',[100,100,900,500])
saveas(gcf,['figures' filesep 'testSimulMetrics_UniAll'],'png')

for simul=1:length(allSig)
    figure;
    subplot(1,3,1+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(aucOld(:,:,simul)),std(aucOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucNew(:,:,simul)),std(aucNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucOldUni(:,:,simul)),std(aucOldUni(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucNewUni(:,:,simul)),std(aucNewUni(:,:,simul)))
    xlabel('nb of sub')
    ylabel('AUC')
    ylim([0 1])
    title(['simul ' nameSim{simul}])
    subplot(1,3,2+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(energyOld(:,:,simul)),std(energyOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyNew(:,:,simul)),std(energyNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyOldUni(:,:,simul)),std(energyOldUni(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyNewUni(:,:,simul)),std(energyNewUni(:,:,simul)))
    xlabel('nb of sub')
    ylabel('Energy')
    ylim([0 1])
    subplot(1,3,3+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(mseOld(:,:,simul)),std(mseOld(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseNew(:,:,simul)),std(mseNew(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseOldUni(:,:,simul)),std(mseOldUni(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseNewUni(:,:,simul)),std(mseNewUni(:,:,simul)))
    xlabel('nb of sub')
    ylabel('MSE')
    ylim([0 2])
    legend('minNorm','newMethod','minNormUni','newMethodUni')
end
set(gcf,'position',[100,100,900,500])
saveas(gcf,['figures' filesep 'testSimulMetrics_Uni'],'png')
