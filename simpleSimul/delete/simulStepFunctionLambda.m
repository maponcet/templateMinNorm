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
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
% dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
noiseLevel = [1 100]; % 1 is same amount of signal/noise, 100 is 100 times more signal than noise
lambda = [5 10 50:50:1000 2000:100:5000 8000 10000]; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

% nbSbjToInclude =[1 2 5 10 20 30 40 50];

totSbj = 1;
totBoot = 20;

% initialise variables
aucAve = zeros(totBoot,length(noiseLevel),length(lambda));
energyAve = zeros(totBoot,length(noiseLevel),length(lambda));
mseAveNorm = zeros(totBoot,length(noiseLevel),length(lambda));

aucWhole = zeros(totBoot,length(noiseLevel),length(lambda));
aucROI = zeros(totBoot,length(noiseLevel),length(lambda));
energyWhole = zeros(totBoot,length(noiseLevel),length(lambda));
energyROI = zeros(totBoot,length(noiseLevel),length(lambda));
mseWholeNorm = zeros(totBoot,length(noiseLevel),length(lambda));
mseROINorm = zeros(totBoot,length(noiseLevel),length(lambda));

numSubs = 20;


for repBoot=1:totBoot
    
% list of random sbj with replacement
listSub = randi(length(dirList),numSubs,1);
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

%% Simulate sources (sourceERP)
% 45ms baseline + 45ms signal 15ms V1 15ms MT 15ms both
x = 0 : pi / 45 : 2 * pi-pi/45; % 360 deg with point every 4 deg
srcAmp = zeros(numROIs,1);
srcERP = zeros(numROIs, length(x) );
% V1 left & right
srcAmp(ac_sources) = 1;
srcERP([ac_sources(1) ac_sources(3)],[46:60 76:90]) = 1;
srcERP([ac_sources(2) ac_sources(4)],[61:90]) = 1;
winERP = 46:90;
% ERP baseline timewindow
timeBase = 1:45;

%%
for nn=1:length(noiseLevel)
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
    [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel(nn),winERP );
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

% figure;plotOnEgi(squeeze(mean(Y_avg(:,:,80),1)));colorbar
% figure;plotOnEgi(squeeze(mean(Y(:,:,80),1)));colorbar
% % % figure;plotOnEgi(squeeze(mean(Y_avgALL(:,:,80),1)));colorbar
% % figure;plotOnEgi(squeeze(Y_avg(1,:,80)));colorbar
% % figure;plotOnEgi(squeeze(Y(1,:,80)));colorbar
% timepoint = 50;
% figure;
% for iSub=1:numSubs
%     subplot(4,5,iSub)
%     plotOnEgi(Y_avg(iSub,:,timepoint)); colorbar
% end
% figure;
% for iSub=1:numSubs
%     subplot(4,5,iSub)
%     plotOnEgi(Y(iSub,:,timepoint)); colorbar
% end



for ll=1:length(lambda)
%% compute minimum norm
regionWhole = zeros(numSubs,numROIs,length(srcERP));
regionROI = zeros(numSubs,numROIs,length(srcERP));
% min_norm on average data: get beta values for each ROI over time
[betaAverage] = minimum_normFix(avMap, squeeze(mean(Y_avg,1)), lambda(ll));


for iSub=1:numSubs
    % regular minimum_norm: on the 20484 indexes per sbj
    [betaWhole] = minimum_normFix(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), lambda(ll));
    yWhole(iSub,:,:) = [fullFwd{iSub}] * betaWhole;

%     [betaWhole, ~, lambdaWhole(iSub,level)] = minimum_norm(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
    %             [betaWhole_SSVEP, ~, lambdaWhole_SSVEP] = minimum_norm(fullFwd{iSub}, squeeze(Y_SSVEPlo(iSub,:,:)), nLambdaRidge);
    %             [betaWhole, ~, lambdaWhole] = minimum_norm(fullFwd{iSub}, squeeze(unstackedY(:,iSub,:)), nLambdaRidge);
    % min_norm on only the ROI indexes per sbj
    [betaROI] = minimum_normFix([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), lambda(ll));
    yROI(iSub,:,:) = [roiFwd{iSub,:}] * betaROI;
%     [betaROI, ~, lambdaROI(iSub,level)] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
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


%% TOPO

% compare with average roiFwd for the participants then plot the average
% betas on this average roiFwd
roiMap = zeros(size(roiFwd{1},1),length(listROIs),numSubs);
for iSub=1:numSubs
    for rr=1:length(listROIs)
        roiMap(:,rr,iSub) = sum([roiFwd{iSub,rr}],2);
    end
end
meanRoiMap = mean(roiMap,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute auc, mse, relative energy using average signal in rois
% do for all the min norm outputs
[aucAve(repBoot,nn,ll), energyAve(repBoot,nn,ll),...
    mseAveNorm(repBoot,nn,ll),] = computeMetrics(betaAverage(:,winERP),srcERP(:,winERP));
[aucWhole(repBoot,nn,ll), energyWhole(repBoot,nn,ll),...
    mseWholeNorm(repBoot,nn,ll)] = computeMetrics(retrieveWhole(:,winERP),srcERP(:,winERP));
[aucROI(repBoot,nn,ll), energyROI(repBoot,nn,ll),...
    mseROINorm(repBoot,nn,ll)] = computeMetrics(retrieveROI(:,winERP),srcERP(:,winERP));



% %% plot BETAs for 1st bootstrap average sbj & ERPs for V1 only
% count = 1;
% figure;set(gcf,'position',[100,100,800,1000])
% for iRoi = 1:2:numROIs
%     % need to normalise the signal
%     subplot(3,3,count);hold on
%     plot(srcERP(iRoi,:) / max(max(abs(srcERP))) )
%     plot(betaAverage(iRoi,:) / max(max(abs(betaAverage))) )
%     plot(retrieveWhole(iRoi,:) / max(max(abs(retrieveWhole))) )
%     plot(retrieveROI(iRoi,:) / max(max(abs(retrieveROI))) )
%     legend('Source','temp','whole','roi','location','best')
%     title(listROIs(iRoi))
%     count=count+1;
% end
% % saveas(gcf,['figures/betasN' num2str(level)],'png')
% 

%%
end

end

end
% save('ERPtestSNR.mat','aucAve','energyAve','mseAveNorm','aucWhole','energyWhole','mseWholeNorm',...
%     'aucROI','energyROI','mseROINorm')



%%% plot metrics
lineCOL={':r',':b',':g','-r','-b','-g'};
figure;
subplot(1,3,1);hold on;
for noise=1:2
plot(lambda,squeeze(mean(aucAve(:,noise,:),1)),lineCOL{1+(noise-1)*3},'LineWidth',2)
plot(lambda,squeeze(mean(aucWhole(:,noise,:),1)),lineCOL{2+(noise-1)*3},'LineWidth',2)
plot(lambda,squeeze(mean(aucROI(:,noise,:),1)),lineCOL{3+(noise-1)*3},'LineWidth',2)
end
xlabel('lambda')
ylabel('AUC')
ylim([0 1])
subplot(1,3,2);hold on;
for noise=1:2
plot(lambda,squeeze(mean(energyAve(:,noise,:),1)),lineCOL{1+(noise-1)*3},'LineWidth',2)
plot(lambda,squeeze(mean(energyWhole(:,noise,:),1)),lineCOL{2+(noise-1)*3},'LineWidth',2)
plot(lambda,squeeze(mean(energyROI(:,noise,:),1)),lineCOL{3+(noise-1)*3},'LineWidth',2)
end
xlabel('lambda')
ylabel('energy')
ylim([0 1])
subplot(1,3,3);hold on;
for noise=1:2
plot(lambda,squeeze(mean(mseAveNorm(:,noise,:),1)),lineCOL{1+(noise-1)*3},'LineWidth',2)
plot(lambda,squeeze(mean(mseWholeNorm(:,noise,:),1)),lineCOL{2+(noise-1)*3},'LineWidth',2)
plot(lambda,squeeze(mean(mseROINorm(:,noise,:),1)),lineCOL{3+(noise-1)*3},'LineWidth',2)
end
xlabel('lambda')
ylabel('mse')
ylim([0 1])
legend('TemplateSNR1','WholeSNR1','ROISNR1','TemplateSNR100','WholeSNR100','ROISNR100')
saveas(gcf,'figures/lambdas','png')

