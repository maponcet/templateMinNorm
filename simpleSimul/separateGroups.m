% simulate EEG for 25 individuals and retrieve the sources for either
% a) templates from the same 25 individuals (same forward)
% b) templates from a different set of 25 individuals
% c) templates used in other simulations (50 individuals)
% source is ERP from one random ROI

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);
numSubs = length(dirList);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 20; % for calculating minimum_norm, hyper param in min norm
totBoot = 60;
aucAve = zeros(totBoot,length(SNRlevel),3);
energyAve = aucAve; mseAveNorm = aucAve;


% ERP & baseline timewindow
timeBase = 1:45;
winERP = 46:90;

for repBoot=1:totBoot
    fprintf('Bootstrap %d\n',repBoot)

    % separateGroups
    randNum = randperm(numSubs);
    gpSig = randNum(1:numSubs/2);
    gpLeft = randNum(numSubs/2+1:end);

    %% LOAD FWD for gpSig & create templates
    fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
    roiMap = zeros(128,numROIs,length(gpLeft));
    for iSub=1:length(gpSig)
        clear fwdMatrix roiInfo
        % fwd file
        load([dataPath dirList(gpSig(iSub)).name])
        fullFwd{iSub} = fwdMatrix;
        % go through each ROI and save the corresponding fwdMesh values
        % corresponding to the indexes of that ROI
        for rr=1:numROIs
            indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            roiMap(:,rr,iSub) = sum(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
            roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
            % to get roiFwd for one sbj= [roiFwd{iSub,:}]
            % save the index for each ROI
            idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
        end
    end
    avMapSameGp = mean(roiMap,3);


    %% create templates for the other 25
    roiMap = zeros(128,numROIs,length(gpLeft));
    for ff=1:length(gpLeft)
        clear fwdMatrix roiInfo
        load([dirList(gpLeft(ff)).folder filesep dirList(gpLeft(ff)).name])
        % go through each ROI and sum all the mesh indexes corresponding to that ROI
        for rr=1:length(listROIs)
            clear indexROI
            indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            roiMap(:,rr,ff) = sum(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
        end
    end
    avMapDiffGp = mean(roiMap,3);



    %% simulate and retrieve signal
    % pick an ROI
    roiNum = randi(numROIs);
    activeROIs = listROIs(roiNum);
    % find the ROI index corresponding to the activeROIs
    ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));
    % create ERP
    srcERP = zeros(numROIs,45*2);
    srcERP(:,winERP) = createSourceERP(numROIs,ac_sources(1));

    for level=1:length(SNRlevel)
        beta = zeros(numROIs,length(srcERP),3);
        noiseLevel = SNRlevel(level);

        %%% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj
        % (using individual fwd model)
        Y = zeros(numSubs/2,size(fullFwd{1},1),length(srcERP));
        Y_avg = zeros(numSubs/2,size(fullFwd{1},1),length(srcERP));

        for iSub=1:numSubs/2
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
        for iSub=1:numSubs/2
            Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
        end

        % retrieve sources
        beta(:,:,1) = minNormFast_lcurve(avMap, squeeze(mean(Y_avg,1)));
        beta(:,:,2) = minNormFast_lcurve(avMapSameGp, squeeze(mean(Y_avg,1)));
        beta(:,:,3)= minNormFast_lcurve(avMapDiffGp, squeeze(mean(Y_avg,1)));
        % compute metrics
        for bb=1:3
            [aucAve(repBoot,level,bb), energyAve(repBoot,level,bb),mseAveNorm(repBoot,level,bb)] = ...
                computeMetrics(beta(:,winERP,bb),srcERP(:,winERP));
        end

    end % noise

end
save('simulOutput/separateGroups.mat','aucAve','energyAve','mseAveNorm')


figure;
for bb=1:3
    subplot(1,3,1);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(aucAve(:,:,bb))),squeeze(std(aucAve(:,:,bb),1)),'LineWidth',2,'CapSize',0)
    xlabel('log10(SNR)');ylim([0 1]);xlim([-1.5 4.5]);ylabel('AUC')
    subplot(1,3,2);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(energyAve(:,:,bb))),squeeze(std(energyAve(:,:,bb),1)),'LineWidth',2,'CapSize',0)
    ylim([0 1]);xlim([-1.5 4.5]);ylabel('Energy');
    subplot(1,3,3);hold on;
    errorbar(log10(SNRlevel),squeeze(mean(mseAveNorm(:,:,bb))),squeeze(std(mseAveNorm(:,:,bb),1)),'LineWidth',2,'CapSize',0)
    ylabel('MSE');ylim([0 1]);xlim([-1.5 4.5]);
end
legend('50','same','diff')
set(gcf,'position',[100 100 1000 300])
saveas(gcf,['figures' filesep 'separateGroups'],'png')
saveas(gcf,['figures' filesep 'separateGroups'],'fig')
print(gcf,['figures' filesep 'separateGroups'],'-depsc')
