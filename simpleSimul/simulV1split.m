clearvars;close all;
% simulate ERP from V1v or V1d L/R and check that it recovers V1 (not other ROI)

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128_allROI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
listROIs = [listROIs {'V1V-L'} {'V1V-R'} {'V1D-L'} {'V1D-R'}];
numROIs = length(listROIs);

% load forward of sbj with distinction between V1v and V1d
nbSub = 0;
for iSub=1:length(dirList)
    clear fwdMatrix roiInfo
    % fwd file
    load([dataPath dirList(iSub).name])
    if ~isempty(find(contains({roiInfo.name},'V1V-L')))
        nbSub = nbSub+1;
        fullFwd{nbSub} = fwdMatrix;
        % go through each ROI and save the corresponding fwdMesh values
        % corresponding to the indexes of that ROI
        for rr=1:numROIs
            indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            roiFwd{nbSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
            idxROIfwd{nbSub,rr} = roiInfo(indexROI).meshIndices;
        end
    end
end

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm
potentialSrc = {'V1V-L','V1V-R','V1D-L','V1D-R'};

totBoot = 30;

for repBoot=1:totBoot
    fprintf('bootstrap %d\n',repBoot)

    %% Simulate sources
    % pick 1 unilateral source in V1
    activeROIs = potentialSrc(randi(4,1));
    % find the ROI index corresponding to the activeROIs
    ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

    % ERP & baseline timewindow
    timeBase = 1:45;
    winERP = 46:90;
    srcERP = zeros(numROIs,45*2); % 45*2 timepoints
    srcERP(:,winERP) = createSourceERP(numROIs,ac_sources(1));

    % list of random sbj with replacement
    listSub = randi(nbSub,nbSub,1);

        for level=1:length(SNRlevel)
            noiseLevel = SNRlevel(level);
            
            %%% Simulate scalp activity (Y)
            % use the generated sources to simulate scalp activity for each sbj
            % (using individual fwd model)
            Y = zeros(nbSub,size(fullFwd{1},1),length(srcERP));
            Y_avg = zeros(nbSub,size(fullFwd{1},1),length(srcERP));
            
            for iSub=1:length(listSub)
                % initialise matrix of source activity
                sourceData = zeros(size(fullFwd{listSub(iSub)},2) , length(srcERP));
                for ss=1:length(ac_sources)
                    sourceData(idxROIfwd{listSub(iSub),ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{listSub(iSub),ac_sources(ss)}),1);
                end
                % multiply fwd (128*20484) with the activated idx over time
                % (sourceData of 20484*90) and obtain Y elec x time
                y_stim = fullFwd{listSub(iSub)} * sourceData;
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
            for iSub=1:nbSub
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
%             [betaLCFUN, betaAverage, betaBest, lambda, lambdaCurv, lambdaBest, ...
%                 lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(Y_avg,1)),srcERP);
            
%             indFwdROI_noise=[roiFwd{iSub,:}];
%             indData_noise=squeeze(Y_avg(iSub,:,:));
%             indFwd_noise=fullFwd{iSub};
%             
%             for iSub=1:numSubs
%                 % regular minimum_norm: on the 20484 indexes per sbj
%                 [betaWhole,lambdaWhole] = minNormFast(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%                 [betaWholeLC, lambdaWholeLC] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
%                 
%                 [betaROI, lambdaROI] = minNormFast([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%                 [betaROILC, lambdaROILC] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
%                 
%                 % beta values are for the indexes, but needs it per ROI
%                 % get the number of indexes per ROI for this subj
%                 rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
%                 % get the range
%                 range = [0 cumsum(rangeROI)]; % cumulative sum of elements
%                 % SUM (not average) the beta values per ROI (=across the indexes)
%                 regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%                 regionROILC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROILC(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%                 
%                 % need to find the indexes for whole brain -> use idxROIfwd
%                 % (no need to get the range)
%                 regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%                 regionWholeLC(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWholeLC(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
%                 
%                 % feed ROI per sbj instead of mesh = "oracle" (best possible recovery)
%                 sbjROI = cell2mat(arrayfun(@(x) sum(fullFwd{iSub}(:,idxROIfwd{iSub,x}),2),1:numROIs,'uni',false));
%                 [betaROIin(iSub,:,:), lambdaGridMinNormROIin] = minNormFast(sbjROI, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
%                 [betaROIinLC(iSub,:,:), lambdaGridMinNormROIinLC] = minNormFast_lcurve(sbjROI, squeeze(Y_avg(iSub,:,:)));
%             end
%             % average across subj
%             retrieveWhole = squeeze(mean(regionWhole,1));
%             retrieveROI = squeeze(mean(regionROI,1));
%             retrieveROIin = squeeze(mean(betaROIin,1));
%             retrieveWholeLC = squeeze(mean(regionWholeLC,1));
%             retrieveROILC = squeeze(mean(regionROILC,1));
%             retrieveROIinLC = squeeze(mean(betaROIinLC,1));
            
            % save simulation
            simulERP(repBoot,level).listROIs = listROIs;
            simulERP(repBoot,level).listSub = listSub;
            simulERP(repBoot,level).winERP = winERP;
            simulERP(repBoot,level).srcERP = srcERP;
            simulERP(repBoot,level).data = Y_avg;
            simulERP(repBoot,level).noise = SNRlevel(level);
            simulERP(repBoot,level).beta = betaAverage;
%             simulERP(totSbj,level).beta(2,:,:) = retrieveWhole;
%             simulERP(totSbj,level).beta(3,:,:) = retrieveROI;
%             simulERP(totSbj,level).beta(4,:,:) = retrieveROIin;
%             simulERP(totSbj,level).beta(5,:,:) = retrieveWholeLC;
%             simulERP(totSbj,level).beta(6,:,:) = retrieveROILC;
%             simulERP(totSbj,level).beta(7,:,:) = retrieveROIinLC;
%             simulERP(totSbj,level).beta(8,:,:) = betaBest;
            
        end % noise
        

end % boot

save('simulOutput/simulV1split.mat','simulERP')




%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
aucAve = zeros(totBoot,length(SNRlevel));energyAve =aucAve ;mseAveNorm=aucAve;
for repBoot = 1:totBoot
    for level=1:length(SNRlevel)
        % srcERP is for subpart of V1 (ventral/dorsal) so change it to V1
        source = zeros(18,length(winERP));
        source(1,:) = sum(simulERP(repBoot,level).srcERP([19 21],winERP)); % left V1
        source(2,:) = sum(simulERP(repBoot,level).srcERP([20 22],winERP)); % right V1
        [aucAve(repBoot,level), energyAve(repBoot,level),mseAveNorm(repBoot,level)] = ...
            computeMetrics(squeeze(simulERP(repBoot,level).beta(:,winERP)),source);
    end
end

figure;hold on
subplot(1,3,1);hold on;
errorbar(log10(SNRlevel),mean(aucAve),std(aucAve),'LineWidth',2,'CapSize',0)
xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
subplot(1,3,2);hold on;
errorbar(log10(SNRlevel),mean(energyAve),std(energyAve),'LineWidth',2,'CapSize',0)
ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
subplot(1,3,3);hold on;
errorbar(log10(SNRlevel),mean(mseAveNorm),std(mseAveNorm),'LineWidth',2,'CapSize',0)
ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
set(gcf,'position',[100 100 800 300])
saveas(gcf,['figures' filesep 'splitV1'],'png')

