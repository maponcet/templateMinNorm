clearvars;close all;

addpath([pwd filesep 'subfunctions' filesep]);

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

load('averageMapEGI32.mat'); % load average map of ROIs (128 elec x 18 ROIs)
mm(1).avMap = avMap;
load('averageMapEGI64.mat'); % load average map of ROIs (128 elec x 18 ROIs)
mm(2).avMap = avMap;
load('averageMapEGI128.mat'); % load average map of ROIs (128 elec x 18 ROIs)
mm(3).avMap = avMap;
load('averageMapEGI256.mat'); % load average map of ROIs (128 elec x 18 ROIs)
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

totBoot = 10;

% initialise variables
aucAve = zeros(totBoot,nbModels,length(SNRlevel));
energyAve = zeros(totBoot,nbModels,length(SNRlevel));
mseAveNorm = zeros(totBoot,nbModels,length(SNRlevel));

aucWhole = zeros(totBoot,nbModels,length(SNRlevel));
aucROI = zeros(totBoot,nbModels,length(SNRlevel));
energyWhole = zeros(totBoot,nbModels,length(SNRlevel));
energyROI = zeros(totBoot,nbModels,length(SNRlevel));
mseWholeNorm = zeros(totBoot,nbModels,length(SNRlevel));
mseROINorm = zeros(totBoot,nbModels,length(SNRlevel));

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
        % min_norm on average data: get beta values for each ROI over time
        % stack all the participants vertically
        stackY = reshape(permute(Y_avg,[2,1,3]),[size(Y_avg,1)*size(Y_avg,2),size(Y_avg,3)]);
        % stack the template for as many participants
        stackAvMap = repmat(mm(sysNb).avMap,numSubs,1);
        [betaAverage, betaMinNorm, lambda, gcvErrorMinNorm, lambdaGridMinNorm] = minimum_norm(stackAvMap, stackY, nLambdaRidge);
             
        for iSub=1:numSubs
            % regular minimum_norm: on the 20484 indexes per sbj
            [betaWhole, betaMinNormWhole, lambdaWhole,...
                gcvErrorMinNormWhole, lambdaGridMinNormWhole] = minimum_norm(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
            % only using the indexes that are in the ROIs
            [betaROI, betaMinNormROI, lambdaROI,...
                gcvErrorMinNormROI, lambdaGridMinNormROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
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
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole,1));
        retrieveROI = squeeze(mean(regionROI,1));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% compute auc, mse, relative energy using average signal in rois
        % do for all the min norm outputs
        [aucAve(repBoot,sysNb,level), energyAve(repBoot,sysNb,level),...
            mseAveNorm(repBoot,sysNb,level),] = computeMetrics(betaAverage(:,winERP),srcERP(:,winERP));
        
        [aucWhole(repBoot,sysNb,level), energyWhole(repBoot,sysNb,level),...
            mseWholeNorm(repBoot,sysNb,level)] = computeMetrics(retrieveWhole(:,winERP),srcERP(:,winERP));
        
        [aucROI(repBoot,sysNb,level), energyROI(repBoot,sysNb,level),...
            mseROINorm(repBoot,sysNb,level)] = computeMetrics(retrieveROI(:,winERP),srcERP(:,winERP));

%         %% Plots for 1st bootstrap
%         if repBoot==1 && numSubs >10
%             %%% plot BETAs
%             count = 1;
%             figure;set(gcf,'position',[100,100,800,1000])            
%             for iRoi = 1:2:numROIs
%                 % need to normalise the signal
%                 subplot(3,3,count);hold on
%                 plot(betaAverage(iRoi,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
%                 plot(betaAverage(iRoi+1,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
%                 tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
%                 ylim([-1 1]);count=count+1;
%             end
%             saveas(gcf,['figures/compareSys' num2str(sysNb)],'png')
%             
%         end
        
    end
    end
    
end

save('compSys.mat','aucAve','energyAve','mseAveNorm','aucROI',...
    'energyROI','mseROINorm','aucWhole','energyWhole','mseWholeNorm')

figure;hold on
for ss=1:sysNb
    subplot(1,3,1);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucAve(:,ss,:))),squeeze(std(aucAve(:,ss,:),1)),'LineWidth',2)
    xlabel('log(SNR)')
    ylim([0 1])
    ylabel('AUC')
    title('template')
    subplot(1,3,2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyAve(:,ss,:))),squeeze(std(energyAve(:,ss,:),1)),'LineWidth',2)
    ylim([0 1])
    ylabel('Energy')
    subplot(1,3,3);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:))),squeeze(std(mseAveNorm(:,ss,:),1)),'LineWidth',2)
    ylabel('MSE')
    ylim([0 1])
end
legend('32','64','124','256','location','best')
set(gcf,'position', [50, 100, 800, 400])
saveas(gcf,['figures/compSysTemplate'],'png')

figure;hold on
for ss=1:sysNb
    subplot(1,3,1);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucROI(:,ss,:))),squeeze(std(aucROI(:,ss,:),1)),'LineWidth',2)
    xlabel('log(SNR)')
    ylim([0 1])
    ylabel('AUC')
    title('ROI')
    subplot(1,3,2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyROI(:,ss,:))),squeeze(std(energyROI(:,ss,:),1)),'LineWidth',2)
    ylim([0 1])
    ylabel('Energy')
    subplot(1,3,3);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseROINorm(:,ss,:))),squeeze(std(mseROINorm(:,ss,:),1)),'LineWidth',2)
    ylabel('MSE')
    ylim([0 1])
end
legend('32','64','124','256','location','best')
set(gcf,'position', [50, 100, 800, 400])
saveas(gcf,['figures/compSysROI'],'png')

figure;hold on
for ss=1:sysNb
    subplot(1,3,1);hold on;
    errorbar(log(SNRlevel),squeeze(mean(aucWhole(:,ss,:))),squeeze(std(aucWhole(:,ss,:),1)),'LineWidth',2)
    xlabel('log(SNR)')
    ylim([0 1])
    ylabel('AUC')
    title('whole')
    subplot(1,3,2);hold on;
    errorbar(log(SNRlevel),squeeze(mean(energyWhole(:,ss,:))),squeeze(std(energyWhole(:,ss,:),1)),'LineWidth',2)
    ylim([0 1])
    ylabel('Energy')
    subplot(1,3,3);hold on;
    errorbar(log(SNRlevel),squeeze(mean(mseWholeNorm(:,ss,:))),squeeze(std(mseWholeNorm(:,ss,:),1)),'LineWidth',2)
    ylabel('MSE')
    ylim([0 1])
end
legend('32','64','124','256','location','best')
set(gcf,'position', [50, 100, 800, 400])
saveas(gcf,['figures/compSysWhole'],'png')


