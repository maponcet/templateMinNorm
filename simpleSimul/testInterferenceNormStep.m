clearvars;close all;
% Test for an interfering source not in visual areas
% simulate ERP from V1+MT 
% retrieve the sources using different methods (template, ROI, whole brain)
% Each simulation uses a different set of sbj and ERP but the same one is
% tested across different levels of noise and number of sbj

% 1-45 = baseline
% 46-90 = V1+MT + interference

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128_allROI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
% sourceL = {'V1-L','MT-L'};
% sourceR = {'V1-R','MT-R'};
sourceL = {'V1-L'};
% sourceR = {'V1-R'};
% sourceExt = {'caudalanteriorcingulate-L','caudalanteriorcingulate-R'};
sourceExt = {'caudalanteriorcingulate-L'};
% simulated signal
% activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
activeROIs = [sourceL]; 
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

nbSbjToInclude =[2 8 20 50];

totBoot = 30;

for repBoot=1:totBoot
    
    %% Simulate sources
    % amplitude (1 to 10) and time function is different for each
    % source but the same for all sbj for a given bootstrap
    % 1-45 = baseline
    % 46-90 = V1-MT
    srcERP = zeros(numROIs,20); % 45*2 timepoints
    srcExtERP = zeros(length(sourceExt),20);
%     srcERP(:,46:90) = createSourceERP(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
    srcERP(ac_sources,:) = 1;
    srcExtERP(1,:) = 1; 
    
%     % ERP & baseline timewindow
%     timeBase = 1:45;
%     winERP = 46:90;
    winERP = 1:20;
        
    for totSbj=1:length(nbSbjToInclude) 
        numSubs = nbSbjToInclude(totSbj);
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
            roiNormVis = vecnorm(fwdMatrix(:,idxROIfwd{iSub,1})'); 
            roiNormExt = vecnorm(fwdMatrix(:,idxExtROIfwd{iSub,1})');
%             figure;plot(roiNormVis);hold on; plot(roiNorExt)
            scalingFactor = roiNormVis / roiNormExt; % mrdivide
            fwdMatrixN = fwdMatrix; fwdMatrixN2 = fwdMatrix; 
            fwdMatrixN(:,idxExtROIfwd{iSub,1}) = scalingFactor * fwdMatrix(:,idxExtROIfwd{iSub,1}); 
            normFwd{iSub} = fwdMatrixN;
            fwdMatrixN2(:,idxExtROIfwd{iSub,1}) = scalingFactor/2 * fwdMatrix(:,idxExtROIfwd{iSub,1}); 
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
                y_stimExt = normFwd{iSub} * sourceDataExt;
                y_stimExt50 = normFwd50{iSub} * sourceDataExt;
                % add gaussian noise on electrodes
                noise = randn(size(y_stimExt,1),size(y_stimExt,2));
                Y(iSub,:,:) = y_stim + (mean(y_stim,2))/sqrt(noiseLevel).*noise ;
                Yext(iSub,:,:) = y_stimExt + (mean(y_stimExt,2))/sqrt(noiseLevel).*noise ;
                Yext50(iSub,:,:) = y_stimExt50 + (mean(y_stimExt50,2))/sqrt(noiseLevel).*noise ;
            
%                 % check SNR
%                 allSig = Y(iSub,:,:); noise1 = (mean(y_stim,2))/sqrt(noiseLevel)*noise;
%                 (rms(allSig(:))/rms(noise1(:)))^2 - 1
%                 allSig = Yext(iSub,:,:); noise1 = (mean(y_stimExt,2))/sqrt(noiseLevel)*noise;
%                 (rms(allSig(:))/rms(noise1(:)))^2 - 1
%                 allSig = Yext50(iSub,:,:); noise1 = (mean(y_stimExt50,2))/sqrt(noiseLevel)*noise;
%                 (rms(allSig(:))/rms(noise1(:)))^2 - 1

            end
            
            %%% Use average reference for centering Y
            %%% that is: substract the average electrode activity at each time point
            % this is done by bsxfun which applies element-wise substraction (the 90
            % averages across electrodes) - useless?
            for iSub=1:numSubs
                Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
                Y_avgExt(iSub,:,:) = bsxfun(@minus,squeeze(Yext(iSub,:,:)), mean(squeeze(Yext(iSub,:,:))));
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
            [betaAverageExt, lambdaExt] = minNormFast_lcurve(avMap, squeeze(mean(Y_avgExt,1)));
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
            simulERP(repBoot,totSbj,level).listROIs = listROIs;
            simulERP(repBoot,totSbj,level).listSub = listSub;
            simulERP(repBoot,totSbj,level).winERP = winERP;
            simulERP(repBoot,totSbj,level).srcERP = srcERP;
            simulERP(repBoot,totSbj,level).data = Y_avg;
            simulERP(repBoot,totSbj,level).dataExt = Y_avgExt50;
            simulERP(repBoot,totSbj,level).noise = SNRlevel(level);
            simulERP(repBoot,totSbj,level).beta(1,:,:) = betaAverage;
            simulERP(repBoot,totSbj,level).beta(2,:,:) = betaAverageExt;
            simulERP(repBoot,totSbj,level).beta(3,:,:) = betaAverageExt50;
%             simulERP(repBoot,totSbj,level).beta(4,:,:) = retrieveWhole;
%             simulERP(repBoot,totSbj,level).beta(5,:,:) = retrieveROI;
%             simulERP(repBoot,totSbj,level).beta(6,:,:) = retrieveROIin;
%             simulERP(repBoot,totSbj,level).beta(5,:,:) = retrieveWholeLC;
%             simulERP(repBoot,totSbj,level).beta(6,:,:) = retrieveROILC;
%             simulERP(repBoot,totSbj,level).beta(7,:,:) = retrieveROIinLC;
            
        end % noise
        
    end % sbj
end % boot

save('simulOutput/simulTestNormExt1srcStep.mat','simulERP','-v7.3')



nbModel = 3;
% initialise variables
aucAve = zeros(size(simulERP,1),size(simulERP,2),size(simulERP,3),nbModel);
energyAve = aucAve;
mseAveNorm = aucAve;

for model=1:nbModel
for repBoot=1:size(simulERP,1)
    for totSbj=1:size(simulERP,2)
        for level=1:size(simulERP,3)            
        [aucAve(repBoot,totSbj,level,model), energyAve(repBoot,totSbj,level,model),mseAveNorm(repBoot,totSbj,level,model)] = ...
            computeMetrics(squeeze(simulERP(repBoot,totSbj,level).beta(model,:,winERP)),simulERP(repBoot,totSbj,level).srcERP(:,winERP));        
        end
    end
end
end

%%% plot metrics
modName = {'template','interference','50interference','whole','ROI','oracle'};
figure;hold on
for model=1:nbModel
for ss=1:size(simulERP,2)
    subplot(3,nbModel,model);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,ss,:,model))),squeeze(std(aucAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    xlabel('log(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
    title(modName(model))
    subplot(3,nbModel,model+nbModel);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,ss,:,model))),squeeze(std(energyAve(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(3,nbModel,model+nbModel*2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,ss,:,model))),squeeze(std(mseAveNorm(:,ss,:,model),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
end
end
legend('N=2','N=8','N=20','N=50')
set(gcf,'position',[100 100 700 700])
saveas(gcf,['figures' filesep 'simulTestNormExt1srcStep'],'png')


% figure;plot(srcERP(1,:)),hold on; plot(srcExtERP(1,:))
% addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
% addpath /Users/marleneponcet/Documents/Git/svndl_code/alesToolbox
% addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
% elecLayout = 'GSN-HydroCel-128.sfp';
% time=57;
% figure;subplot(2,4,1);plotTopo(mean(Y_avg(:,:,time)),elecLayout);colorbar
% subplot(2,4,2);plotTopo(mean(Y_avgExt(:,:,time)),elecLayout);colorbar
% subplot(2,4,3);plotTopo(mean(Y_avgExt50(:,:,time)),elecLayout);colorbar
% subplot(2,4,4);imagesc([betaAverage(:,time) betaAverageExt(:,time) betaAverageExt50(:,time)]);colorbar
% time=48;
% subplot(2,4,5);plotTopo(mean(Y_avg(:,:,time)),elecLayout);colorbar
% subplot(2,4,6);plotTopo(mean(Y_avgExt(:,:,time)),elecLayout);colorbar
% subplot(2,4,7);plotTopo(mean(Y_avgExt50(:,:,time)),elecLayout);colorbar
% subplot(2,4,8);imagesc([betaAverage(:,time) betaAverageExt(:,time) betaAverageExt50(:,time)]);colorbar
% 
% figure;plot(squeeze(mean(Y_avgExt(:,1,:))))
% hold on; plot(squeeze(mean(Y_avgExt50(:,1,:))))