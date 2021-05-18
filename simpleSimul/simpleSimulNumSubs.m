clearvars;close all;
% simulation1: V1+MT
% simulation2: V2v+V4
% Each source has a random amplitude and waveform but is the same across
% participants (also change for different simulations)

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
allSig = {{'V1-R','V1-L','MT-R','MT-L'}};%{'V2V-R','V4-R','V2V-L','V4-L'}}; % all areas simulated (could try unilateral with total 9 areas?)
nbSbjToInclude =[1 2 4 8 16 25 32 50];
totBoot = 20; % nb of bootstrap

% initialise variables
aucMinNorm = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
aucGpMinNorm = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
aucMinNormInd = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyMinNorm = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyGpNorm = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
energyMinNormInd = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseMinNorm = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseGpNorm = zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
mseMinNormInd= zeros(totBoot,length(nbSbjToInclude),length(allSig)); 
    
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
            x = 0 : pi / 45 : 2 * pi-pi/45;
            sourceValOverTime = zeros(numROIs, length(x) );
            sourceValOverTimeERP = zeros(numROIs, length(x) );
            for k = 1 : length(ac_sources)
                % random amplitude between 1 and 10
                sourceAmplitude = 1 + randi(90)/10;
                % create SSVEP waveform with random phase delay and random 
                % coef for each source
                y = rand * cos( 4 * x - rand * pi ) + ...
                    rand * cos( 8 * x - rand * pi ) + ...
                    rand * cos( 12 * x - rand * pi ) + ...
                    rand * cos( 16 * x - rand * pi );
                % create ERP waveform (use cos envelop instead?)
                yERP = y.*gaussmf(x,[pi/5 pi]);      
                % multiply waveform by the source amplitude (ROIxtime)
                sourceValOverTime( ac_sources(k) , : ) = y * sourceAmplitude;
                sourceValOverTimeERP( ac_sources(k) , : ) = yERP * sourceAmplitude;                
            end

            % Simul EEG signal from sourceValOverTime for each sbj (different
            % ROI map)
            for iSub=1:numSubs
                y_stim = roiMap(:,:,iSub) * sourceValOverTime;
                y_stimERP = roiMap(:,:,iSub) * sourceValOverTimeERP;
                % add noise
                [noisy_data,noise_var] = add_noise_with_SNR( y_stim , SNRlevel );
                Y(128*(iSub-1)+1:128*iSub,:) = y_stim + noisy_data;          
                [noisy_data,noise_var] = add_noise_with_SNR( y_stimERP , SNRlevel );
                Yerp(128*(iSub-1)+1:128*iSub,:) = y_stimERP + noisy_data; 
            end
            
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
            [betaReg, ~, lambdaReg] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
            [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(roiStacked, Ylo, nLambdaRidge);
            [betaRegERP, ~, lambdaRegERP] = minimum_norm(avMap.activity, grandMeanDataERP, nLambdaRidge);
            [betaMinNormERP, ~, lambdaMinNormERP] = minimum_norm(roiStacked, YloERP, nLambdaRidge);
            % betaMinNorm is numSubs*numROIs x time so need to extract betaMinNorm of a single sbj
            % for each ROI separately then average across sbj
            regionTmp = arrayfun(@(x) betaMinNorm(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
            regionMinNorm =reshape(cell2mat(regionTmp),numROIs,size(Ylo,2),numSubs);
            regionMinNormAvg = mean(regionMinNorm,3);
            regionTmpERP = arrayfun(@(x) betaMinNormERP(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
            regionMinNormERP =reshape(cell2mat(regionTmpERP),numROIs,size(Ylo,2),numSubs);
            regionMinNormAvgERP = mean(regionMinNormERP,3);            
            
            %% compute auc, mse, relative energy using average signal in rois
            %%%%% AUC
            tmp = zeros(numROIs,1);
            tmp(ac_sources,:) = 1;
            % instead of using 1 and 0 in the matrix, use the real amount of simulated
            % activity? Nooo rocArea only uses 1 and 0s
            % AUC computed without normalising: it wouldn't change the results
            for nT = 1:size(Ylo,2)
                aucGpMinNormT(nT) = rocArea( abs(betaReg(:,nT)) , tmp );
                aucMinNormT(nT) = rocArea( abs(regionMinNormAvg(:,nT)) , tmp );
                aucGpMinNormT_ERP(nT) = rocArea( abs(betaRegERP(:,nT)) , tmp );
                aucMinNormT_ERP(nT) = rocArea( abs(regionMinNormAvgERP(:,nT)) , tmp );                
            end
            aucGpMinNorm(repBoot,totSbj,simulation) = mean(aucGpMinNormT);
            aucMinNorm(repBoot,totSbj,simulation) = mean(aucMinNormT); % auc calculated from the average activity (beta)
            aucGpMinNorm_ERP(repBoot,totSbj,simulation) = mean(aucGpMinNormT_ERP);
            aucMinNorm_ERP(repBoot,totSbj,simulation) = mean(aucMinNormT_ERP); 
            
            %%%%% relative energy
            for nT = 1:size(Ylo,2)
                norm_betaReg(:,nT) = betaReg(:,nT) / max( abs(betaReg(:,nT)) ); % normalise estimated sources
                relEnergy(nT) = sum( abs( norm_betaReg(ac_sources,nT) ) ) / sum( abs(norm_betaReg(:,nT)) );
                norm_MinNorm(:,nT) = regionMinNormAvg(:,nT) / max( abs(regionMinNormAvg(:,nT)) ); % normalise estimated sources
                relEnergyMinNorm(nT) = sum( abs( norm_MinNorm(ac_sources,nT) ) ) / sum( abs(norm_MinNorm(:,nT)) );
                norm_betaRegERP(:,nT) = betaRegERP(:,nT) / max( abs(betaRegERP(:,nT)) ); % normalise estimated sources
                relEnergyERP(nT) = sum( abs( norm_betaRegERP(ac_sources,nT) ) ) / sum( abs(norm_betaRegERP(:,nT)) );
                norm_MinNormERP(:,nT) = regionMinNormAvgERP(:,nT) / max( abs(regionMinNormAvgERP(:,nT)) ); % normalise estimated sources
                relEnergyMinNormERP(nT) = sum( abs( norm_MinNormERP(ac_sources,nT) ) ) / sum( abs(norm_MinNormERP(:,nT)) );
           end
            energyGpNorm(repBoot,totSbj,simulation) = mean(relEnergy);
            energyMinNorm(repBoot,totSbj,simulation) = mean(relEnergyMinNorm);
            energyGpNorm_ERP(repBoot,totSbj,simulation) = mean(relEnergyERP);
            energyMinNorm_ERP(repBoot,totSbj,simulation) = mean(relEnergyMinNormERP);
            
            %%%%% mse
            % paper: norm(sourceValOverTime(:,nT)-betaReg(:,nT),'fro').^2 / numel(betaReg(:,nT))
            for nT = 1:size(Ylo,2)
                normSourceTime = sourceValOverTime(:,nT) / max(abs(sourceValOverTime(:,nT))); % normalise simulated source
                mseBetaReg(nT) = sum( (normSourceTime(ac_sources) - norm_betaReg(ac_sources,nT)).^2 ) / sum( (normSourceTime(ac_sources)).^2 );
                mseMin(nT) = sum( (normSourceTime(ac_sources) - norm_MinNorm(ac_sources,nT)).^2 ) / sum( (normSourceTime(ac_sources)).^2 );
                normSourceTimeERP = sourceValOverTimeERP(:,nT) / max(abs(sourceValOverTimeERP(:,nT))); % normalise simulated source
                mseBetaRegERP(nT) = sum( (normSourceTimeERP(ac_sources) - norm_betaRegERP(ac_sources,nT)).^2 ) / sum( (normSourceTimeERP(ac_sources)).^2 );
                mseMinERP(nT) = sum( (normSourceTimeERP(ac_sources) - norm_MinNormERP(ac_sources,nT)).^2 ) / sum( (normSourceTimeERP(ac_sources)).^2 );
            end
            mseGpNorm(repBoot,totSbj,simulation) = mean(mseBetaReg);
            mseMinNorm(repBoot,totSbj,simulation) = mean(mseMin);
            mseGpNorm_ERP(repBoot,totSbj,simulation) = mean(mseBetaRegERP);
            mseMinNorm_ERP(repBoot,totSbj,simulation) = mean(mseMinERP);
            
            %% plot source amplitude overtime per ROI
            % for 1)true signal 2)new template minnorm 3)standard minnorm
            if numSubs == 50 && repBoot == 1
                figure;rr=1;
                set(gcf,'position',[100,100,800,1000])
                for iRoi = 1:2:size(betaReg,1)
                    subplot(9,3,1+3*(rr-1));hold on;
                    plot(sourceValOverTime(iRoi,:)','m-','linewidth',2)
                    plot(sourceValOverTime(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(sourceValOverTime)) max(max(sourceValOverTime))])
                    title(listROIs(iRoi))                    
                    subplot(9,3,2+3*(rr-1));hold on;
                    plot(betaReg(iRoi,:)','m-','linewidth',2)
                    plot(betaReg(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(betaReg)) max(max(betaReg))])
                    subplot(9,3,3+3*(rr-1));hold on;
                    plot(regionMinNormAvg(iRoi,:)','m-','linewidth',2)
                    plot(regionMinNormAvg(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(regionMinNormAvg)) max(max(regionMinNormAvg))])
                    rr=rr+1;
                end
                saveas(gcf,['figures' filesep 'testSimul' num2str(simulation)],'png')
                figure;rr=1;
                set(gcf,'position',[100,100,800,1000])
                for iRoi = 1:2:size(betaReg,1)
                    subplot(9,3,1+3*(rr-1));hold on;
                    plot(sourceValOverTimeERP(iRoi,:)','m-','linewidth',2)
                    plot(sourceValOverTimeERP(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(sourceValOverTimeERP)) max(max(sourceValOverTimeERP))])
                    title(listROIs(iRoi))
                    subplot(9,3,2+3*(rr-1));hold on;
                    plot(betaRegERP(iRoi,:)','m-','linewidth',2)
                    plot(betaRegERP(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(betaRegERP)) max(max(betaRegERP))])
                    subplot(9,3,3+3*(rr-1));hold on;
                    plot(regionMinNormAvgERP(iRoi,:)','m-','linewidth',2)
                    plot(regionMinNormAvgERP(iRoi+1,:)','b-','linewidth',2)
                    ylim([min(min(regionMinNormAvgERP)) max(max(regionMinNormAvgERP))])
                    rr=rr+1;
                end
                saveas(gcf,['figures' filesep 'testSimulERP' num2str(simulation)],'png')
            end
            
            %% difference in computing metrics
            % I could compute metrics per sbj (only for regular MinNorm) or on
            % average... completely different results...
            clear aucMinNormTSub norm_MinNormI normSourceTimeI relEnergyMinNormI
            for ss=1:numSubs
                for nT = 1:size(Ylo,2)
                    aucMinNormTSub(nT,ss) = rocArea( abs(regionMinNorm(:,nT,ss)) , tmp );
                end
            end
            aucMinNormInd(repBoot,totSbj,simulation) = mean(mean(aucMinNormTSub)); % auc calculated for each sbj then average
            for ss=1:numSubs
                for nT = 1:size(Ylo,2)
                    norm_MinNormI(:,nT,ss) = regionMinNorm(:,nT,ss) / max( abs(regionMinNorm(:,nT,ss)) ); % normalise estimated sources
                    relEnergyMinNormI(nT,ss) = sum( abs( norm_MinNormI(ac_sources,nT,ss) ) ) / sum( abs(norm_MinNormI(:,nT,ss)) );
                end
            end
            energyMinNormInd(repBoot,totSbj,simulation) = mean(mean(relEnergyMinNormI));
            % not sure how to normalise - here by sbj but could have been for the
            % entire group...
            for ss=1:numSubs
                for nT = 1:size(Ylo,2)
                    normSourceTimeI = sourceValOverTime(ac_sources,nT) / max(abs(sourceValOverTime(ac_sources,nT))); % normalise simulated source
                    n_indMSE(nT,ss) = mean(sum( (normSourceTimeI - norm_MinNormI(ac_sources,nT,ss)).^2 ) / sum( (normSourceTimeI).^2 ));
                end
            end
            mseMinNormInd(repBoot,totSbj,simulation) = mean(mean(n_indMSE));
           
        end % end boostrap
        
    end % go through numSub
    
end % different activated sources


nameSim = {'V1-MT','V2V-V4'};
for simul=1:length(allSig)
    figure;
    subplot(1,3,1+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(aucMinNorm(:,:,simul)),std(aucMinNorm(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucGpMinNorm(:,:,simul)),std(aucGpMinNorm(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucMinNormInd(:,:,simul)),std(aucMinNormInd(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucMinNorm_ERP(:,:,simul)),std(aucMinNorm_ERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(aucGpMinNorm_ERP(:,:,simul)),std(aucGpMinNorm_ERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('AUC')
    ylim([0 1])
    title(['simul ' nameSim{simul}])
    subplot(1,3,2+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(energyMinNorm(:,:,simul)),std(energyMinNorm(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyGpNorm(:,:,simul)),std(energyGpNorm(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyMinNormInd(:,:,simul)),std(energyMinNormInd(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyMinNorm_ERP(:,:,simul)),std(energyMinNorm_ERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(energyGpNorm_ERP(:,:,simul)),std(energyGpNorm_ERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('Energy')
    ylim([0 1])
    subplot(1,3,3+(simul-1)*3);hold on;
    errorbar(nbSbjToInclude,mean(mseMinNorm(:,:,simul)),std(mseMinNorm(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseGpNorm(:,:,simul)),std(mseGpNorm(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseMinNormInd(:,:,simul)),std(mseMinNormInd(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseMinNorm_ERP(:,:,simul)),std(mseMinNorm_ERP(:,:,simul)))
    errorbar(nbSbjToInclude,mean(mseGpNorm_ERP(:,:,simul)),std(mseGpNorm_ERP(:,:,simul)))
    xlabel('nb of sub')
    ylabel('MSE')
    ylim([0 2])
    legend('minNorm','newMethod','minNormInd','minNormERP','newMethodERP')
end
set(gcf,'position',[100,100,900,500])
saveas(gcf,['figures' filesep 'testSimulMetrics'],'png')

