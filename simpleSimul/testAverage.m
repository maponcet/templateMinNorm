clearvars;close all;
% Compare min norm on average EEG vs, single ind then average betas

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm

nbSbjToInclude =[2 8 20 50];

totBoot = 30;


aucAve = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyAve = aucAve;mseAveNorm = aucAve;
aucInd = aucAve; energyInd = aucAve;mseIndNorm = aucAve;

for repBoot=1:totBoot
    
    %% Simulate sources
    roiActive = randsample(1:2:numROIs,2);
    sourceL = listROIs(roiActive);
    sourceR = listROIs(roiActive+1);
    activeROIs = [sourceL,sourceR];
    % find the ROI index corresponding to the activeROIs
    ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));
    
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
                % multiply fwd (128*20484) with the activated idx over time
                % (sourceData of 20484*90) and obtain Y elec x time
                y_stim = fullFwd{iSub} * sourceData;
                % add noise
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
            end
            
            %% compute minimum norm
            [betaAverage, lambdaAv] = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
            % per ind 
            clear beta
            for iSub=1:numSubs
                [beta(iSub,:,:), lambda(iSub)] = minNormFast_lcurve(avMap, squeeze(Y_avg(iSub,:,:)));
            end
            betaInd = squeeze(mean(beta));
            % stack ind
            clear Ystack avMapStack betaStack
            % stack sbj vertically
            test = permute(Y_avg,[2 1 3]);
            Ystack = reshape(test,[size(Y_avg,2)*numSubs,size(Y_avg,3)]);
            avMapStack = repmat(avMap,numSubs,1);
            [betaStack, lambdaStack] = minNormFast_lcurve(avMapStack, Ystack);            
            
            % figure;imagesc(betaInd)
            % figure;imagesc(betaAverage)
            
            [aucAve(repBoot,totSbj,level), energyAve(repBoot,totSbj,level), mseAveNorm(repBoot,totSbj,level)] = ...
                computeMetrics(squeeze(betaAverage(:,winERP)),srcERP(:,winERP));
            [aucInd(repBoot,totSbj,level), energyInd(repBoot,totSbj,level), mseIndNorm(repBoot,totSbj,level)] = ...
                computeMetrics(squeeze(betaInd(:,winERP)),srcERP(:,winERP));
            [aucSta(repBoot,totSbj,level), energySta(repBoot,totSbj,level), mseStaNorm(repBoot,totSbj,level)] = ...
                computeMetrics(squeeze(betaStack(:,winERP)),srcERP(:,winERP));            
        end
    end
end

figure;
for ss=length(nbSbjToInclude)
subplot(1,3,1);hold on;
errorbar(log(SNRlevel),squeeze(mean(aucAve(:,ss,:))),squeeze(std(aucAve(:,ss,:),1)),'LineWidth',2)
errorbar(log(SNRlevel),squeeze(mean(aucInd(:,ss,:))),squeeze(std(aucInd(:,ss,:),1)),'LineWidth',2)
errorbar(log(SNRlevel),squeeze(mean(aucSta(:,ss,:))),squeeze(std(aucSta(:,ss,:),1)),'LineWidth',2)
ylabel('AUC');xlabel('SNR');
subplot(1,3,2);hold on;
errorbar(log(SNRlevel),squeeze(mean(energyAve(:,ss,:))),squeeze(std(energyAve(:,ss,:),1)),'LineWidth',2)
errorbar(log(SNRlevel),squeeze(mean(energyInd(:,ss,:))),squeeze(std(energyInd(:,ss,:),1)),'LineWidth',2)
errorbar(log(SNRlevel),squeeze(mean(energySta(:,ss,:))),squeeze(std(energySta(:,ss,:),1)),'LineWidth',2)
ylabel('Energy');xlabel('SNR');
subplot(1,3,3);hold on;
errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:))),squeeze(std(mseAveNorm(:,ss,:),1)),'LineWidth',2)
errorbar(log(SNRlevel),squeeze(mean(mseIndNorm(:,ss,:))),squeeze(std(mseIndNorm(:,ss,:),1)),'LineWidth',2)
errorbar(log(SNRlevel),squeeze(mean(mseStaNorm(:,ss,:))),squeeze(std(mseStaNorm(:,ss,:),1)),'LineWidth',2)
ylabel('MSE');xlabel('SNR');
end
legend('average','ind then average','stacked')

