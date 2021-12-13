% look at the regularisation parameter for each sbj using gcv or l_curve

clearvars;close all

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
% SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
SNRlevel = [0.1 10000];
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L'};
sourceR = {'V1-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

numSubs = 10; listSub = 1:10;
fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
% fwd file
for iSub=1:numSubs
    load([dataPath dirList(listSub(iSub)).name])
    fullFwd{iSub} = fwdMatrix;
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
        % to get roiFwd for one sbj= [roiFwd{iSub,:}]
        % save the index for each ROI
        idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
    end
end

%% Simulate sources
% 45ms baseline + 45ms signal V1
x = 0 : pi / 45 : 2 * pi-pi/45; % 360 deg with point every 4 deg
srcAmp = zeros(numROIs,1);
srcERP = zeros(numROIs, length(x) );
% V1 left & right
srcAmp(ac_sources) = 1;
srcERP([ac_sources(1) ac_sources(2)],[46:90]) = 1;
% srcERP([ac_sources(2) ac_sources(4)],[61:90]) = 1;
winERP = 46:90;
% ERP baseline timewindow
timeBase = 1:45;

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
        for ss=1:length(ac_sources)
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
        Y_noise(iSub,:,:) = noisy_data;
    end
    
    %%% Use average reference for centering Y
    for iSub=1:numSubs
        Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
    end
    
    %% compute minimum norm
    for iSub=1:numSubs
        %         [betaROI, betaMinNormROI, lambdaROI, lambdaGridROI] = minNorm_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
        [beta_W, lambda_W] = minNormFast_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
        %           [betaROI, lambdaROI] = minNormFast_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
        %         % regular minimum_norm: on the 20484 indexes per sbj
        % %         [betaWhole, betaMinNormWhole, lambdaWhole,...
        % %             gcvErrorMinNormWhole, lambdaGridMinNormWhole] = minimum_norm(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
        %         [betaROI, betaMinNormROI, lambdaROI,...
        %             gcvErrorMinNormROI, lambdaGridMinNormROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), nLambdaRidge);
        %
        % beta values are for the indexes, but I want it per ROI
        % get the number of indexes per ROI for this subj
        rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
        % get the range
        range = [0 cumsum(rangeROI)]; % cumulative sum of elements
        % sum the beta values per ROI (=across the indexes)
        %         regionROI(level,iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
        %         regionWhole(level,iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
        
        lambda_ROI(level,iSub)=lambdaROI;
        lambdaGrid_ROI(level,iSub,:,:) = lambdaGridROI;
        
        %         gcvGrid_ROI(level,iSub,:,:) = gcvErrorMinNormROI;
        % %         lambda_Whole(level,iSub)=lambdaWhole;
        % %         lambdaGrid_Whole(level,iSub,:,:) = lambdaGridMinNormWhole;
        % %         gcvGrid_Whole(level,iSub,:,:) = gcvErrorMinNormWhole;
        
        
    end
    
end

figure;
for iSub=1:50
    subplot(5,10,iSub);hold on;
    for level=1%:length(SNRlevel)
        plot(squeeze(lambdaGrid_Whole(level,iSub,:)),squeeze(gcvGrid_Whole(level,iSub,:)))
    end
    %     plot([lambda(level,iSub) lambda(level,iSub)],ylim,'r:','linewidth',2)
    xlabel('lambda')
    ylabel('gcv error')
end
% lambda does not change with noise

figure;
for iSub=1:50
    subplot(5,10,iSub);hold on;
    plot(squeeze(regionROI(2,iSub,:,90)))
    xlabel('ROI')
    ylabel('beta')
end
minNorm_lcurve( stackAvMap, stackY );

[betaROI, betaMinNormROI, lambdaROI, lambdaGridROI] = minNorm_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
% Optimization of the regularization parameter according to the gcv error
[u,s,v] = csvd(fullFwd{2});
%[lambda, gcvErrorMinNorm, lambdaGridMinNorm] = gcv(u,s,b,'Tikh', nLambda);
[lambda,rho,eta,lambdaGridMinNorm] = l_curve(u,s,squeeze(Y_avg(iSub,:,90))','Tikh');







%%%%%% look at outputs for different lambdas
clearvars;close all;

addpath([pwd filesep 'subfunctions' filesep]);
addpath([pwd filesep 'reguTime' filesep]);
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
% SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
SNRlevel = 10000;
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','V4-L','MT-L'};
sourceR = {'V1-R','V4-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

numSubs = 1; listSub = 1:50;
fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
% fwd file
for iSub=1:numSubs
    load([dataPath dirList(listSub(5)).name])
    fullFwd{iSub} = fwdMatrix;
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
        % to get roiFwd for one sbj= [roiFwd{iSub,:}]
        % save the index for each ROI
        idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
    end
end
%% Simulate sources
% 45ms baseline + 45ms signal V1
x = 0 : pi / 45 : 2 * pi-pi/45; % 360 deg with point every 4 deg
srcAmp = zeros(numROIs,1);
srcERP = zeros(numROIs, length(x) );
% V1 left & right
srcAmp(ac_sources) = 1;
srcERP([ac_sources(1) ac_sources(4)],[46:55 76:90]) = 1;
srcERP([ac_sources(2) ac_sources(5)],[56:65 76:90]) = 1;
srcERP([ac_sources(3) ac_sources(6)],[66:90]) = 1;
winERP = 46:90;
% ERP baseline timewindow
timeBase = 1:45;

noiseLevel = SNRlevel;
%%% Simulate scalp activity (Y)
% use the generated sources to simulate scalp activity for each sbj
% (using individual fwd model)
Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
Y_noise = Y;
Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));

% initialise matrix of source activity
sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
sourceNoise = sourceData; % no signal, used to compute SNR
for ss=1:length(ac_sources)
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
Y_noise(iSub,:,:) = noisy_data;


%%% Use average reference for centering Y
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end

%% compute minimum norm
[beta_gcv, betaGrid_gcv, lambda_gcv, gcvError_gcv,lambdaGrid_gcv] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), 30);
[betaROI, betaMinNormROI, lambdaROI, lambdaGridROI] = minNorm_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
range = [0 cumsum(rangeROI)]; % cumulative sum of elements

for lambda=1:length(lambdaGridROI)
    clear regionROI 
    currentBeta = betaMinNormROI{lambda};
    regionROI = cell2mat(arrayfun(@(x) sum(currentBeta(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    [aucAve(lambda), energyAve(lambda), mseAve(lambda)] = computeMetrics(regionROI(:,winERP),srcERP(:,winERP)); 
end

for lambda=1:length(gcvError_gcv)
    clear regionGCV
    currentBeta_gcv = betaGrid_gcv{lambda};
    regionGCV = cell2mat(arrayfun(@(x) sum(currentBeta_gcv(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    [aucAve_gcv(lambda), energyAve_gcv(lambda), mseAve_gcv(lambda)] = computeMetrics(regionGCV(:,winERP),srcERP(:,winERP));   
end

figure;
subplot(1,3,1)
semilogx(lambdaGridROI,aucAve,'lineWidth',2);hold on;
semilogx(lambdaGrid_gcv,aucAve_gcv,'lineWidth',2)
semilogx([lambdaROI lambdaROI],ylim,'b:','linewidth',2)
semilogx([lambda_gcv lambda_gcv],ylim,'r:','linewidth',2)
xlabel('log lambda');ylabel('AUC')
subplot(1,3,2);semilogx(lambdaGridROI,energyAve,'lineWidth',2);hold on;
semilogx(lambdaGrid_gcv,energyAve_gcv,'lineWidth',2)
semilogx([lambdaROI lambdaROI],ylim,'b:','linewidth',2)
semilogx([lambda_gcv lambda_gcv],ylim,'r:','linewidth',2)
xlabel('log lambda');ylabel('energy')
legend('lcurve','gcv','location','best')
subplot(1,3,3);semilogx(lambdaGridROI,mseAve,'lineWidth',2);hold on;
semilogx(lambdaGrid_gcv,mseAve_gcv,'lineWidth',2)
semilogx([lambdaROI lambdaROI],ylim,'b:','linewidth',2)
semilogx([lambda_gcv lambda_gcv],ylim,'r:','linewidth',2)
xlabel('log lambda');ylabel('MSE')

[beta_w, betaMinNorm_w, lambda_w, lambdaGrid_w] = minNorm_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
[beta_wgcv, betaMinNorm_wgcv, lambda_wgcv, gcvError_wgcv,lambdaGrid_wgcv] = minimum_norm(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)),30);
for iLambda = 1:length(lambdaGrid_w)
    mse_w(iLambda) = sum(sum((betaMinNorm_w{iLambda}-sourceData).^2));
end
for iLambda = 1:length(lambdaGrid_wgcv)
    mse_wgcv(iLambda) = sum(sum((betaMinNorm_wgcv{iLambda}-sourceData).^2));
end
figure;
loglog(lambdaGrid_w,mse_w,'linewidth',2)
hold on;
loglog(lambdaGrid_wgcv,mse_wgcv,'linewidth',2)
loglog([lambda_w lambda_w],ylim,'b--','linewidth',2)
loglog([lambda_wgcv lambda_wgcv],ylim,'r--','linewidth',2)
ylabel('log MSE to ground truth')
xlabel('log Lambda')

[beta_fix] = minimum_normFix(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), 0);
regionFix = cell2mat(arrayfun(@(x) sum(beta_fix(idxROIfwd{iSub,x},:)),1:numROIs,'uni',false)');
[beta_fixROI] = minimum_normFix([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), 0);
regionROIfix = cell2mat(arrayfun(@(x) sum(beta_fixROI(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
[aucAve_fix, energyAve_fix, mseAve_fix] = computeMetrics(regionROIfix(:,winERP),srcERP(:,winERP));
mse_fix = sum(sum((beta_fix-sourceData).^2));

results = [[lambdaGrid_gcv(end) lambdaGridROI(end) 0]
[energyAve_gcv(end) energyAve(end) energyAve_fix]
[mseAve_gcv(end) mseAve(end) mseAve_fix]
[NaN mse_w(end) mse_fix]]';








%%%%%% compare gcv curvature vs lcfun
clearvars;close all;

addpath(genpath([pwd filesep 'subfunctions']))
addpath(['previousFunMinNorm' filesep])
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 10000];
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

%% Simulate sources
x = 0 : pi / 45 : 2 * pi-pi/45; % 360 deg with point every 4 deg
srcAmp = zeros(numROIs,1);
srcERP = zeros(numROIs,45*4); % 45*4 timepoints
srcERP(:,46:90) = createSourceERP(numROIs,ac_sources(1),ac_sources(3));
srcERP(:,91:135) = createSourceERP(numROIs,ac_sources(2),ac_sources(4));
srcERP(:,136:180) = createSourceERP(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
% ERP baseline timewindow
timeBase = 1:45;winERP = 46:90;

numSubs = 1; listSub = randi(50);
fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
% fwd file
for iSub=1:numSubs
    load([dataPath dirList(listSub(iSub)).name])
    fullFwd{iSub} = fwdMatrix;
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
    Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
    Y_noise = Y;
    Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
    
    % initialise matrix of source activity
    sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
    sourceNoise = sourceData; % no signal, used to compute SNR
    for ss=1:length(ac_sources)
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
    Y_noise(iSub,:,:) = noisy_data;
    
    %%% Use average reference for centering Y
    for iSub=1:numSubs
        Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
    end
    
    %% compute minimum norm
%     [beta_gcvROI, betaGrid_gcvROI, lambda_gcvROI, gcvError_gcvROI,lambdaGrid_gcvROI] = minimum_norm([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)), 30);
%     [betaROI, betaMinNormROI, lambdaROI, lambdaGridROI, lambdaCurveROI] = minNorm_lcurve([roiFwd{iSub,:}], squeeze(Y_avg(iSub,:,:)));
     clear mse_gcv 
    [beta_gcv, betaGrid_gcv, lambda_gcv, gcvError_gcv,lambdaGrid_gcv] = minimum_norm(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)), 10);
    [beta, betaMinNorm, lambda, lambdaGrid, lambdaCurve] = minNorm_lcurve(fullFwd{iSub}, squeeze(Y_avg(iSub,:,:)));
    for iLambda = 1:length(lambdaGrid_gcv)
        mse_gcv(iLambda) = sum(sum((betaGrid_gcv{iLambda}-sourceData).^2));
    end
%     for iLambda = 1:length(lambdaGridROI)
%         mse(iLambda) = sum(sum((betaMinNormROI{iLambda}-sourceData).^2));
%     end
   
    figure;
    semilogx(lambdaGrid_gcv,mse_gcv,'linewidth',2)
    hold on;
    plot([lambda lambda],ylim,'b--','linewidth',2)
    plot([lambdaCurve lambdaCurve],ylim,'g--','linewidth',2)
    plot([lambda_gcv lambda_gcv],ylim,'r--','linewidth',2)
    ylabel('MSE to ground truth')
    xlabel('log Lambda')
    legend('','lcfun','curvature','gcv')

end
