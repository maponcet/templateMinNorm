clearvars;close all;
% test 10 sbj, template 50, template 40 different, template 10 same sbj

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
% dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

totBoot = 10;
numSubs = 10; % nb of sbj included in the simulated signal

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
    fprintf('bootstrap %d\n',repBoot)
    % list of random sbj without replacement
    listSub = randperm(length(dirList),numSubs);
    otherSbj = setdiff(1:length(dirList),listSub);
    % since everything is fixed should sample without replacement -> RANDPERM
    % (otherwise same sbj means that it's twice the same)

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
% fprintf('noise %d\n',level)
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

%%% Use average reference for centering Y
%%% that is: substract the average electrode activity at each time point
% this is done by bsxfun which applies element-wise substraction (the 90
% averages across electrodes) - Useless
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end


%% compute minimum norm using 3 different templates
% stack all the participants vertically
stackY = reshape(permute(Y_avg,[2,1,3]),[size(Y_avg,1)*size(Y_avg,2),size(Y_avg,3)]);

% compute the template for the same sbj 
for iSub=1:numSubs
    for rr=1:length(listROIs)
        roiMapOverlap(:,rr,iSub) = sum(roiFwd{iSub,rr},2); 
    end
end
avMapOverlap = mean(roiMapOverlap,3);

% no overlap
for iSub=1:length(otherSbj)
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
% for rr=1:length(listROIs)
%     figure;
%     subplot(1,3,1);plotOnEgi(avMap(:,rr));colorbar
%     subplot(1,3,2);plotOnEgi(avMapOverlap(:,rr));colorbar
%     subplot(1,3,3);plotOnEgi(avMapNoOverlap(:,rr));colorbar
% end

% compute min norm
% 50 sbj templates
% stack the template for as many participants
stackAvMap = repmat(avMap,numSubs,1);
[betaAverage] = minimum_norm(stackAvMap, stackY, nLambdaRidge);
stackOverlap = repmat(avMapOverlap,numSubs,1);
[betaAverageOverlap] = minimum_norm(stackOverlap, stackY, nLambdaRidge);
stackNoOverlap = repmat(avMapNoOverlap,numSubs,1);
[betaAverageNoOverlap] = minimum_norm(stackNoOverlap, stackY, nLambdaRidge);

[testA] = minimum_norm(avMapNoOverlap, squeeze(mean(Y_avg)), nLambdaRidge);

% for iSub=1:numSubs
%     [testB(:,:,iSub)] = minimum_norm(avMapNoOverlap, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
% end
% tt = mean(testB,3);
% figure;imagesc(tt);
% figure;imagesc(betaAverage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute auc, mse, relative energy using average signal in rois
% do for all the min norm outputs
[aucAve(repBoot,level), energyAve(repBoot,level), mseAve(repBoot,level),] = computeMetrics(betaAverage(:,winERP),srcERP(:,winERP));
[aucAveO(repBoot,level), energyAveO(repBoot,level), mseAveO(repBoot,level),] = computeMetrics(betaAverageOverlap(:,winERP),srcERP(:,winERP));
[aucAveNO(repBoot,level), energyAveNO(repBoot,level), mseAveNO(repBoot,level),] = computeMetrics(betaAverageNoOverlap(:,winERP),srcERP(:,winERP));



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

legend('T50','T10overlap','T40noOverlap')
saveas(gcf,['figures' filesep 'templateOverlapNorm2'],'png')
savefig(['figures' filesep 'templateOverlapNorm2.fig'])


