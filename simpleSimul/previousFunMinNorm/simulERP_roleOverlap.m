clearvars;close all;
% test the role of sbj overlap between template and sbj
% only use the template method
% 3 conditions: 25test+25template full overlap, 25+25 no overlap,
% 25test+50template (template that is used in other analysis)
% test for several levels of SNR

addpath([pwd filesep 'subfunctions' filesep]);
addpath([pwd filesep 'reguTime' filesep]);
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
% dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

totBoot = 10;
numSubs = 25; % nb of sbj included in the simulated signal

% initialise variables
aucAve = zeros(totBoot,length(SNRlevel));
energyAve = zeros(totBoot,length(SNRlevel));
mseAve = zeros(totBoot,length(SNRlevel));
aucAveO = zeros(totBoot,length(SNRlevel));
energyAveO = zeros(totBoot,length(SNRlevel));
mseAveO = zeros(totBoot,length(SNRlevel));
aucAveNO = zeros(totBoot,length(SNRlevel));
energyAveNO = zeros(totBoot,length(SNRlevel));
mseAveNO = zeros(totBoot,length(SNRlevel));


for repBoot=1:totBoot
%     fprintf('bootstrap %d\n',repBoot)
    % list of random sbj without replacement -> RANDPERM
    listSub = randperm(length(dirList),numSubs);
    otherSbj = setdiff(1:length(dirList),listSub);

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
fprintf('bootstrap %d noise %d\n',repBoot,level)
%%
noiseLevel = SNRlevel(level);

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
 
% averages across electrodes) - Useless
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end


%% compute minimum norm using 3 different templates
% stack all the participants vertically
stackY = reshape(permute(Y_avg,[2,1,3]),[size(Y_avg,1)*size(Y_avg,2),size(Y_avg,3)]);

% compute the template for the 25 sbj overlap 
for iSub=1:numSubs
    for rr=1:length(listROIs)
        roiMapOverlap(:,rr,iSub) = sum(roiFwd{iSub,rr},2); 
    end
end
avMapOverlap = mean(roiMapOverlap,3);

% 25 no overlap
for iSub=1:numSubs
    clear fwdMatrix roiInfo
    load([dataPath dirList(otherSbj(iSub)).name])
    for rr=1:numROIs
        clear idxROI
        idxROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiMapNoOverlap(:,rr,iSub) = sum(fwdMatrix(:,roiInfo(idxROI).meshIndices),2);
    end
end
avMapNoOverlap = mean(roiMapNoOverlap,3);

% % Compare the ROI templates = they are identical but different scales!
% figure;
% for rr=1:length(listROIs)
%     subplot(3,6,rr);plotOnEgi(avMapOverlap(:,rr));colorbar
% end
% figure;
% for rr=1:length(listROIs)
%     subplot(3,6,rr);plotOnEgi(avMapNoOverlap(:,rr));colorbar
% end
% figure;subplot(3,1,1);imagesc(avMapOverlap);colorbar;title('roi overlap')
% subplot(3,1,2);imagesc(avMapNoOverlap);colorbar;title('roi no overlap')
% subplot(3,1,3);imagesc(avMapOverlap-avMapNoOverlap);colorbar;title('difference')

% compute min norm
% 50 sbj templates
% stack the template for as many participants
stackAvMap = repmat(avMap,numSubs,1);
[betaAverage, lambda] = minNormFast_lcurve(stackAvMap, stackY);
stackOverlap = repmat(avMapOverlap,numSubs,1);
[betaAverageOverlap, lambdaO] = minNormFast_lcurve(stackOverlap, stackY);
stackNO = repmat(avMapNoOverlap,numSubs,1);
[betaAverageNO, lambdaNO] = minNormFast_lcurve(stackNO, stackY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute auc, mse, relative energy using average signal in rois
% do for all the min norm outputs
[aucAve(repBoot,level), energyAve(repBoot,level), mseAve(repBoot,level),] = computeMetrics(betaAverage(:,winERP),srcERP(:,winERP));
[aucAveO(repBoot,level), energyAveO(repBoot,level), mseAveO(repBoot,level),] = computeMetrics(betaAverageOverlap(:,winERP),srcERP(:,winERP));
[aucAveNO(repBoot,level), energyAveNO(repBoot,level), mseAveNO(repBoot,level),] = computeMetrics(betaAverageNO(:,winERP),srcERP(:,winERP));



%%%%% test for scaled betas
[betaAverage_scaled] = minNormFast_lcurve_scaled(stackAvMap, stackY);
[betaAverageOverlap_scaled] = minNormFast_lcurve_scaled(stackOverlap, stackY);
[betaAverageNO_scaled] = minNormFast_lcurve_scaled(stackNO, stackY);

[aucAveScaled(repBoot,level), energyAveScaled(repBoot,level), mseAveScaled(repBoot,level),] = computeMetrics(betaAverage_scaled(:,winERP),srcERP(:,winERP));
[aucAveScaledO(repBoot,level), energyAveScaledO(repBoot,level), mseAveScaledO(repBoot,level),] = computeMetrics(betaAverageOverlap_scaled(:,winERP),srcERP(:,winERP));
[aucAveScaledNO(repBoot,level), energyAveScaledNO(repBoot,level), mseAveScaledNO(repBoot,level),] = computeMetrics(betaAverageNO_scaled(:,winERP),srcERP(:,winERP));



% do a normalisation for each ROI -> unit norming
% so that the betas represent microVolts (instead of microVolts/area size
% as it is now)
% unit norming is: all electrodes are squared and summed. These values are
% then divided so that the total of the electrodes for each ROI (power) is
% equal to 1
regParam = sqrt(sum(avMap.^2,1));
avMapNorm = bsxfun(@rdivide,avMap,regParam);
stackAvMapNorm = repmat(avMapNorm,numSubs,1);
[betaNorm] = minNormFast_lcurve(stackAvMapNorm, stackY);

regParam = sqrt(sum(avMapOverlap.^2,1));
avMapNormO = bsxfun(@rdivide,avMapOverlap,regParam);
stackAvMapNormO = repmat(avMapNormO,numSubs,1);
[betaNormOverlap] = minNormFast_lcurve(stackAvMapNormO, stackY);

regParam = sqrt(sum(avMapNoOverlap.^2,1));
avMapNormNO = bsxfun(@rdivide,avMapNoOverlap,regParam);
stackAvMapNormNO = repmat(avMapNormNO,numSubs,1);
[betaNormNoOverlap] = minNormFast_lcurve(stackAvMapNormNO, stackY);

% % Compare the ROI templates = they are identical
% for rr=1:length(listROIs)
%     figure;
%     subplot(1,3,1);plotOnEgi(avMapNorm(:,rr));colorbar
%     subplot(1,3,2);plotOnEgi(avMapNormO(:,rr));colorbar
%     subplot(1,3,3);plotOnEgi(avMapNormNO(:,rr));colorbar
% end

[aucAveNorm(repBoot,level), energyAveNorm(repBoot,level), mseAveNorm(repBoot,level),] = computeMetrics(betaNorm(:,winERP),srcERP(:,winERP));
[aucAveNormO(repBoot,level), energyAveNormO(repBoot,level), mseAveNormO(repBoot,level),] = computeMetrics(betaNormOverlap(:,winERP),srcERP(:,winERP));
[aucAveNormNO(repBoot,level), energyAveNormNO(repBoot,level), mseAveNormNO(repBoot,level),] = computeMetrics(betaNormNoOverlap(:,winERP),srcERP(:,winERP));



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % normalise ?????????
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % normalise templates?
% for iSub=1:numSubs
%     roiMapNO_norm(:,:,iSub) = roiMapNoOverlap(:,:,iSub) ./ (max(max(abs(roiMapNoOverlap(:,:,iSub)))));
% end
% avMapNO_norm = mean(roiMapNO_norm,3);
% 
% % % Compare the ROI templates 
% % for rr=1
% %     figure;
% %     subplot(2,2,1);plotOnEgi(avMap(:,rr));colorbar
% %     subplot(2,2,2);plotOnEgi(avMapOverlap(:,rr));colorbar
% %     subplot(2,2,3);plotOnEgi(avMapNoOverlap(:,rr));colorbar
% %     subplot(2,2,4);plotOnEgi(avMapNO_norm(:,rr));colorbar
% % end
% 
% stackNoOverlap_norm = repmat(avMapNO_norm,numSubs,1);
% % [betaNO_norm] = minimum_norm(stackNoOverlap_norm, stackY, nLambdaRidge);
% % 
% % [aucNO_norm(repBoot,level), energyNO_norm(repBoot,level), mseNO_norm(repBoot,level)] = computeMetrics(betaNO_norm(:,winERP),srcERP(:,winERP));
% 
% % normalise Y?
% for iSub=1:numSubs
%     Y_avgNorm(iSub,:,:) = Y_avg(iSub,:,:) ./ max(max(abs(Y_avg(iSub,:,:)),[],3));
% end
% % stack all the participants vertically
% stackY_norm = reshape(permute(Y_avgNorm,[2,1,3]),[size(Y_avgNorm,1)*size(Y_avgNorm,2),size(Y_avgNorm,3)]);
% [betaNO_normY] = minimum_norm(stackNoOverlap_norm, stackY_norm, nLambdaRidge);
% % figure; imagesc(betaNO_normY)
% % figure; imagesc(betaNO_normY2)
% % figure; imagesc(betaNorm)
% [aa(repBoot,level), ee(repBoot,level), mm(repBoot,level)] = computeMetrics(betaNO_normY(:,winERP),srcERP(:,winERP));


end

end


%%% plot metrics
figure;
subplot(1,3,1);hold on;
errorbar(log(SNRlevel),mean(aucAve),std(aucAve),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveO),std(aucAveO),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveNO),std(aucAveNO),'LineWidth',2)
xlabel('log(SNR)')
ylim([0 1])
ylabel('AUC')

subplot(1,3,2);hold on;
errorbar(log(SNRlevel),mean(energyAve),std(energyAve),'LineWidth',2)
errorbar(log(SNRlevel),mean(energyAveO),std(energyAveO),'LineWidth',2)
errorbar(log(SNRlevel),mean(energyAveNO),std(energyAveNO),'LineWidth',2)
ylim([0 1])
ylabel('Energy')

subplot(1,3,3);hold on;
errorbar(log(SNRlevel),mean(mseAve),std(mseAve),'LineWidth',2)
errorbar(log(SNRlevel),mean(mseAveO),std(mseAveO),'LineWidth',2)
errorbar(log(SNRlevel),mean(mseAveNO),std(mseAveNO),'LineWidth',2)
ylabel('MSE')
ylim([0 1])

legend('T50','T25overlap','T25noOverlap')
saveas(gcf,['figures' filesep 'templateOverlap'],'png')
savefig(['figures' filesep 'templateOverlap.fig'])


figure;
subplot(1,3,1);hold on;
errorbar(log(SNRlevel),mean(aucAveNorm),std(aucAveNorm),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveNormO),std(aucAveNormO),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveNormNO),std(aucAveNormNO),'LineWidth',2)
xlabel('log(SNR)')
ylim([0 1])
ylabel('AUC')

subplot(1,3,2);hold on;
errorbar(log(SNRlevel),mean(energyAveNorm),std(energyAveNorm),'LineWidth',2)
errorbar(log(SNRlevel),mean(energyAveNormO),std(energyAveNormO),'LineWidth',2)
errorbar(log(SNRlevel),mean(energyAveNormNO),std(energyAveNormNO),'LineWidth',2)
ylim([0 1])
ylabel('Energy')

subplot(1,3,3);hold on;
errorbar(log(SNRlevel),mean(mseAveNorm),std(mseAveNorm),'LineWidth',2)
errorbar(log(SNRlevel),mean(mseAveNormO),std(mseAveNormO),'LineWidth',2)
errorbar(log(SNRlevel),mean(mseAveNormNO),std(mseAveNormNO),'LineWidth',2)
ylabel('MSE')
ylim([0 1])

legend('T50','T25overlap','T25noOverlap')
saveas(gcf,['figures' filesep 'templateOverlapNorm'],'png')
savefig(['figures' filesep 'templateOverlapNorm.fig'])


figure;
subplot(1,3,1);hold on;
errorbar(log(SNRlevel),mean(aucAveScaled),std(aucAveScaled),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveScaledO),std(aucAveScaledO),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveScaledNO),std(aucAveScaledNO),'LineWidth',2)
xlabel('log(SNR)')
ylim([0 1])
ylabel('AUC')

subplot(1,3,2);hold on;
errorbar(log(SNRlevel),mean(energyAveScaled),std(energyAveScaled),'LineWidth',2)
errorbar(log(SNRlevel),mean(energyAveScaledO),std(energyAveScaledO),'LineWidth',2)
errorbar(log(SNRlevel),mean(energyAveScaledNO),std(energyAveScaledNO),'LineWidth',2)
ylim([0 1])
ylabel('Energy')

subplot(1,3,3);hold on;
errorbar(log(SNRlevel),mean(mseAveScaled),std(mseAveScaled),'LineWidth',2)
errorbar(log(SNRlevel),mean(mseAveScaledO),std(mseAveScaledO),'LineWidth',2)
errorbar(log(SNRlevel),mean(mseAveScaledNO),std(mseAveScaledNO),'LineWidth',2)
ylabel('MSE')
ylim([0 1])

legend('T50','T25overlap','T25noOverlap')
saveas(gcf,['figures' filesep 'templateOverlapScaled'],'png')
savefig(['figures' filesep 'templateOverlapScaled.fig'])
