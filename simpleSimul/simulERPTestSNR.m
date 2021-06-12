clearvars;close all;
% simulation: V1+MT
% for a given simulation, amplitude and time function is different for each
% source but the same across participants (with different fwd models)
% simulation consistent with retinotopy: L&R sources are the same
% (= assumes full field stimulation)
% simulation on MESH
% minimum_norm: done on a)whole brain, b)only 18 ROIs, c)average50ROIs
% TEST ERP SNR

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 10 50 100 500 1000]; % noise level 10% = if the signal is 10 then the noise is 10*10, for 0.5 S=50/100
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

% nbSbjToInclude =[1 2 5 10 20 30 40 50];
nbSbjToInclude =[3 20 50];
totBoot = 5; % nb of bootstrap

% initialise variables
aucAve = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyAve = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseAveNorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));

aucWhole = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
aucROI = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyWhole = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
energyROI = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseWholeNorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));
mseROINorm = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));

snrRawTime = zeros(totBoot,length(nbSbjToInclude),length(SNRlevel));

for totSbj=1:length(nbSbjToInclude)
    numSubs = nbSbjToInclude(totSbj);
for level=1:length(SNRlevel)
        noiseLevel = SNRlevel(level);   
    for repBoot =1:totBoot
        fprintf('sbj%d noise %d bootstrap %d \n',numSubs,noiseLevel,repBoot);
        
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
        [srcAmp, srcSSVEP, srcERP,winERP] = createSourceROI(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
        % ERP baseline timewindow
        timeBase = setdiff(1:size(srcERP,2),winERP);
        
        %% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj 
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
        Y_noise = Y; 
        
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
            % add noise
            [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
            % to keep the same SNR for the 2 Y, need to compute noise for 
            % the 2 Y separately as it is based on the variance of the signal
            Y(iSub,:,:) = y_stim + noisy_data;
%             Y_SSVEPprev(iSub,:,:) = y_stimSSVEP + noisy_data;
            Y_noise(iSub,:,:) = noisy_data;
%             Y_noiseSSVEP(iSub,:,:) = noisy_dataSSVEP;
%             Ypure(iSub,:,:) = y_stim;
%             Y2(iSub,:,:) = y_stim + noisy_dataSSVEP;
        end 
                
        %% Use average reference for centering Y 
        for iSub=1:numSubs
            Y(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
%             Y2(iSub,:,:) = bsxfun(@minus,squeeze(Y2(iSub,:,:)), mean(squeeze(Y2(iSub,:,:))));
        end
 
        % check SNR
        ySNR_time=zeros(1,numSubs);
        for iSub=1:numSubs
            ySNR_time(iSub) = (rms(rms(Y(iSub,:,winERP)))/rms(rms(Y(iSub,:,timeBase)))) ^2 -1;
        end       
        snrRawTime(repBoot,totSbj,level) = mean(ySNR_time); 
        
        
        %%% in previous code dim reduction per sbj -> results pretty bad
        %%% and nb of sbj did not matter much!
%         %%%% test dim reduction
%         for iSub=1:numSubs
%             [u1, s1, v1] = svd(squeeze(Y(iSub,:,:)));
%             YloIND(iSub,:,:) = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
%         end
%         % stack sbj vertically
%         test = permute(Y,[2 1 3]);
%         test2 = reshape(test,[128*10,90]);
%         n = numel(test2);
%         [u1, s1, v1] = svd(test2);
%         Ylo10 = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
%         unstackedData = reshape(Ylo10,128,10,size(Ylo10,2));
%         grandMeanData = squeeze(mean(unstackedData,2));
%         figure;plotContourOnScalp(squeeze(mean(YloIND(1:5,:,45),1)),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
%         figure;plotContourOnScalp(squeeze(mean(YloIND(1:10,:,45),1)),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
%         figure;plotContourOnScalp(squeeze(grandMeanData(:,45)),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
%         testH = permute(Y(1:5,:,:),[2 1 3]);
%         testH = reshape(testH,[128*5,90]);
%         n = numel(testH);
%         [u1, s1, v1] = svd(testH);
%         YH = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
%         unstackedDataH = reshape(YH,128,5,size(YH,2));
%         grandMeanDataH = squeeze(mean(unstackedDataH,2));
%         figure;plotContourOnScalp(squeeze(grandMeanDataH(:,45)),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')

%         %%  reduce dimension data (Ylo)
%         % =PCA denoised version of Y (denoised by truncation of the SVD)
%         %%%%%% per sbj or alltogether?
%         for iSub=1:numSubs
%             [u1, s1, v1] = svd(squeeze(Y(iSub,:,:)));
%             Ylo(iSub,:,:) = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
%         end
%         % compute SNR after reduction... without -1 because i only have
%         % signal now?
%         snrLo = zeros(1,numSubs);snrLoTime = zeros(1,numSubs);
%         for iSub=1:numSubs
%             snrLo(iSub) = (rms(rms(Ylo(iSub,:,:)))/rms(rms(Y_noise(iSub,:,:)))) ^2 ;
%             snrLoTime(iSub) = (rms(rms(Ylo(iSub,:,winERP)))/rms(rms(Ylo(iSub,:,timeBase)))) ^2 ;            
%         end
%         snrPCA(repBoot,totSbj,level) = mean(snrLo);
%         snrPCAtime(repBoot,totSbj,level) = mean(snrLoTime);
        

        
%        % stack sbj vertically
%        tmpY = permute(Y,[2 1 3]);
%        tmpY2 = reshape(tmpY,[size(fullFwd{1},1)*numSubs,length(srcERP)]);
%        [u1, s1, v1] = svd(tmpY2);
%        YloStack = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
%        tmpYs = permute(Y_SSVEP,[2 1 3]);
%        tmpY2s = reshape(tmpYs,[size(fullFwd{1},1)*numSubs,length(srcERP)]);
%        [u1, s1, v1] = svd(tmpY2s);
%        Y_SSVEPloStack = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
       
        %% compute minimum norm
        regionWhole = zeros(numSubs,numROIs,length(srcERP));
        regionROI = zeros(numSubs,numROIs,length(srcERP));        
        % min_norm on average data: get beta values for each ROI over time 
        [betaAverage, ~, lambdaAverage] = minimum_norm(avMap, squeeze(mean(Y,1)), nLambdaRidge);
%         unstackedY = reshape(YloStack,size(fullFwd{1},1),numSubs,size(YloStack,2));
%         [betaAverageStack, ~, lambdaAverageStack] = minimum_norm(avMap, squeeze(mean(unstackedY,2)), nLambdaRidge);
%         unstackedSSVEP = reshape(Y_SSVEPloStack,size(fullFwd{1},1),numSubs,size(Y_SSVEPloStack,2));
%         [betaAverageSSVEP, ~, lambdaAverageSSVEP] = minimum_norm(avMap, squeeze(mean(unstackedSSVEP,2)), nLambdaRidge);
        for iSub=1:numSubs
            % regular minimum_norm: on the 20484 indexes per sbj
            [betaWhole, ~, lambdaWhole] = minimum_norm(fullFwd{iSub}, squeeze(Y(iSub,:,:)), nLambdaRidge);
%             [betaWhole_SSVEP, ~, lambdaWhole_SSVEP] = minimum_norm(fullFwd{iSub}, squeeze(Y_SSVEPlo(iSub,:,:)), nLambdaRidge);
%             [betaWhole, ~, lambdaWhole] = minimum_norm(fullFwd{iSub}, squeeze(unstackedY(:,iSub,:)), nLambdaRidge);
            % min_norm on only the ROI indexes per sbj
            [betaROI, ~, lambdaROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Y(iSub,:,:)), nLambdaRidge);
%             [betaROI_SSVEP, ~, lambdaROI_SSVEP] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_SSVEPlo(iSub,:,:)), nLambdaRidge);
%             [betaROIStack, ~, lambdaROIStack] = minimum_norm([roiFwd{iSub,:}], squeeze(unstackedY(:,iSub,:)), nLambdaRidge);
            
            % beta values are for the indexes, but I want it per ROI
            % get the number of indexes per ROI for this subj
            rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
            % get the range
            range = [0 cumsum(rangeROI)]; % cumulative sum of elements
            % average the beta values per ROI (=across the indexes)
            regionROI(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); 
%             regionROI_SSVEP(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROI_SSVEP(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); 
%             regionROI_Stack(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaROIStack(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); 
            % need to find the indexes for whole brain
            regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
%             regionWhole_SSVEP(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole_SSVEP(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
            % create brain resp for each sbj
            % YhatWhole(iSub,:,:) = fullFwd{iSub}*betaWhole(iSub,:,:);
        end
        % average across subj
        retrieveWhole = squeeze(mean(regionWhole,1));
        retrieveROI = squeeze(mean(regionROI,1));
%         retrieveROI_SSVEP = squeeze(mean(regionROI_SSVEP,1));
%         retrieveWhole_SSVEP = squeeze(mean(regionWhole_SSVEP,1));
%         retrieveROI_Stack = squeeze(mean(regionROI_Stack,1));
        
        
        
        %% compute auc, mse, relative energy using average signal in rois
        % do for all the min norm outputs
        [aucAve(repBoot,totSbj,level), energyAve(repBoot,totSbj,level),...
            mseAveNorm(repBoot,totSbj,level),] = computeMetrics(betaAverage(:,winERP),ac_sources,srcERP(:,winERP));
        [aucWhole(repBoot,totSbj,level), energyWhole(repBoot,totSbj,level),...
            mseWholeNorm(repBoot,totSbj,level)] = computeMetrics(retrieveWhole(:,winERP),ac_sources,srcERP(:,winERP));
        [aucROI(repBoot,totSbj,level), energyROI(repBoot,totSbj,level),...
            mseROINorm(repBoot,totSbj,level)] = computeMetrics(retrieveROI(:,winERP),ac_sources,srcERP(:,winERP));
        
        
        %% plot BETAs for 1st bootstrap average sbj & ERPs for V1 only
        if repBoot==1
            count = 1;
            figure;set(gcf,'position',[100,100,800,1000])
            for iRoi = 1:2:numROIs
                % need to normalise the signal
                subplot(3,3,count);hold on
                plot(srcERP(iRoi,:) / max(max(abs(srcERP))) )
                plot(betaAverage(iRoi,:) / max(max(abs(betaAverage))) )
                plot(retrieveWhole(iRoi,:) / max(max(abs(retrieveWhole))) )
                plot(retrieveROI(iRoi,:) / max(max(abs(retrieveROI))) )
                title(['SNR ' num2str(snrRawTime(repBoot,totSbj,level),'%.1f') ])
                legend('Source','temp','whole','roi')
                title(listROIs(iRoi))
                count=count+1;
            end
        end
      
      
    end % end boostrap
end % SNR


end % go through nb of sub to include

save('ERPtestSNR.mat','aucAve','energyAve','mseAveNorm','aucWhole','energyWhole','mseWholeNorm',...
    'aucROI','energyROI','mseROINorm')

%%% plot metrics
lineCOL={':r',':b',':g','--r','--b','--g','-r','-b','-g'};
figure;
subplot(1,3,1);hold on;
for ss=1:length(nbSbjToInclude)
errorbar(SNRlevel,squeeze(mean(aucAve(:,ss,:))),squeeze(std(aucAve(:,ss,:))),lineCOL{1+(ss-1)*3})
errorbar(SNRlevel,squeeze(mean(aucWhole(:,ss,:))),squeeze(std(aucWhole(:,ss,:))),lineCOL{2+(ss-1)*3})
errorbar(SNRlevel,squeeze(mean(aucROI(:,ss,:))),squeeze(std(aucROI(:,ss,:))),lineCOL{3+(ss-1)*3})
end
xlabel('SNR (1sbj)')
ylabel('AUC')
ylim([0 1])
subplot(1,3,2);hold on;
for ss=1:length(nbSbjToInclude)
errorbar(SNRlevel,squeeze(mean(energyAve(:,ss,:))),squeeze(std(energyAve(:,ss,:))),lineCOL{1+(ss-1)*3})
errorbar(SNRlevel,squeeze(mean(energyWhole(:,ss,:))),squeeze(std(energyWhole(:,ss,:))),lineCOL{2+(ss-1)*3})
errorbar(SNRlevel,squeeze(mean(energyROI(:,ss,:))),squeeze(std(energyROI(:,ss,:))),lineCOL{3+(ss-1)*3})
end
xlabel('SNR (1sbj)')
ylabel('Energy')
ylim([0 1])
subplot(1,3,3);hold on;
for ss=1:length(nbSbjToInclude)
errorbar(SNRlevel,squeeze(mean(mseAveNorm(:,ss,:))),squeeze(std(mseAveNorm(:,ss,:))),lineCOL{1+(ss-1)*3})
errorbar(SNRlevel,squeeze(mean(mseWholeNorm(:,ss,:))),squeeze(std(mseWholeNorm(:,ss,:))),lineCOL{2+(ss-1)*3})
errorbar(SNRlevel,squeeze(mean(mseROINorm(:,ss,:))),squeeze(std(mseROINorm(:,ss,:))),lineCOL{3+(ss-1)*3})
end
xlabel('SNR (1sbj)')
ylabel('MSEnorm')
ylim([0 1])
legend('3New','3Whole','3ROI','20New','20Whole','20ROI','50New','50Whole','50ROI')
set(gcf,'position',[100,100,1200,400])
saveas(gcf,['figures' filesep 'testERPsnr'],'png')
