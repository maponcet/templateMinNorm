clearvars;close all

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
% SNRlevel = [0.1 1 10 200 10000]; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
SNRlevel = [0.1  10000];
% SNRlevel = 10;
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

%% Simulate sources
% % 45ms baseline + 45ms signal V1 
% x = 0 : pi / 45 : 2 * pi-pi/45; % 360 deg with point every 4 deg
% srcAmp = zeros(numROIs,1);
% srcERP = zeros(numROIs, length(x) );
% % V1 left & right
% srcAmp(ac_sources) = 1;
% srcERP([ac_sources(1) ac_sources(3)],[46:60 76:90]) = 1;
% srcERP([ac_sources(2) ac_sources(4)],[61:75 76:90]) = 1;
% % srcERP([ac_sources(2) ac_sources(4)],[61:90]) = 1;
% winERP = 46:90;
% % ERP baseline timewindow
% timeBase = 1:45;

numSubs = 30;
stackAvMap = repmat(avMap,numSubs,1);

srcERP = zeros(numROIs,45*4); % 45*4 timepoints
srcERP(:,46:90) = createSourceERP(numROIs,ac_sources(1),ac_sources(3));
srcERP(:,91:135) = createSourceERP(numROIs,ac_sources(2),ac_sources(4));
srcERP(:,136:180) = createSourceERP(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));

% ERP & baseline timewindow
timeBase = 1:45;
winERP = 46:180; 


listSub = randi(length(dirList),numSubs,1);
listSub = randperm(50,30); 
listSub = repmat(randperm(50,6),5,1);
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


clear betaBest_Lcurve beta_Lcurve lamda_Lcurve lambdaGrid_Lcurve
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
    % stack all the participants vertically
    stackY = reshape(permute(Y_avg,[2,1,3]),[size(Y_avg,1)*size(Y_avg,2),size(Y_avg,3)]);
    [betaBest_Lcurve, beta_Lcurve, lamda_Lcurve, ...
        lambdaGrid_Lcurve,lambdaCurv] = minNorm_lcurve(  stackAvMap, stackY);
    clear mseTruth
    for iLambda = 1:length(beta_Lcurve)
        mseTruth(iLambda) = sum(sum((beta_Lcurve{iLambda}-srcERP).^2));
    end
    figure;
    semilogy(lambdaGrid_Lcurve,mseTruth,'linewidth',2)
    hold on
    semilogy([lamda_Lcurve lamda_Lcurve],ylim,'g--','linewidth',2)
    semilogy([lambdaCurv lambdaCurv],ylim,'r--','linewidth',2)
    xlabel('lambda');ylabel('log MSE')
    legend('','lcfun','curvature')
    
end

save('sameSbjBadCurvatureBadLcfun','stackAvMap', 'stackY')

level=1;
for iLambda = 1:length(beta_Lcurve)
    mseTruth(iLambda) = sum(sum((beta_Lcurve{level,iLambda}-srcERP).^2));
    [aucAve(iLambda), energyAve(iLambda), mseAve(iLambda)] = computeMetrics(beta_Lcurve{level,iLambda}(:,winERP),srcERP(:,winERP)); 
end
figure;
semilogy(lambdaGrid_Lcurve(level,:),mseTruth,'linewidth',2)
hold on
semilogy([lamda_Lcurve(level) lamda_Lcurve(level)],ylim,'k--','linewidth',2)
xlabel('lambda');ylabel('log MSE')
% saveas(gcf,['figures' filesep 'mseTemplateERP_1'],'png')

% figure;
% subplot(1,3,1)
% semilogx(lambdaGrid_Lcurve(level,:),aucAve,'lineWidth',2);hold on;
% semilogx([lamda_Lcurve(level) lamda_Lcurve(level)],ylim,'b:','linewidth',2)
% xlabel('log lambda');ylabel('AUC')
% subplot(1,3,2);semilogx(lambdaGrid_Lcurve(level,:),energyAve,'lineWidth',2);hold on;
% semilogx([lamda_Lcurve(level) lamda_Lcurve(level)],ylim,'b:','linewidth',2)
% xlabel('log lambda');ylabel('energy')
% subplot(1,3,3);semilogx(lambdaGrid_Lcurve(level,:),mseAve,'lineWidth',2);hold on;
% semilogx([lamda_Lcurve(level) lamda_Lcurve(level)],ylim,'b:','linewidth',2)
% xlabel('log lambda');ylabel('MSE')
% % saveas(gcf,['figures' filesep 'metricsTemplateERPlog_1'],'png')

figure;
subplot(1,3,1)
plot(lambdaGrid_Lcurve(level,:),aucAve,'lineWidth',2);hold on;
plot([lamda_Lcurve(level) lamda_Lcurve(level)],ylim,'b:','linewidth',2)
xlabel('lambda');ylabel('AUC')
subplot(1,3,2);plot(lambdaGrid_Lcurve(level,:),energyAve,'lineWidth',2);hold on;
plot([lamda_Lcurve(level) lamda_Lcurve(level)],ylim,'b:','linewidth',2)
xlabel('lambda');ylabel('energy')
subplot(1,3,3);plot(lambdaGrid_Lcurve(level,:),mseAve,'lineWidth',2);hold on;
plot([lamda_Lcurve(level) lamda_Lcurve(level)],ylim,'b:','linewidth',2)
xlabel('lambda');ylabel('MSE')
% saveas(gcf,['figures' filesep 'metricsTemplateERP_1'],'png')

figure;
subplot(2,1,1);imagesc(squeeze(srcERP));title('source ERP');xlabel('time');ylabel('ROI')
subplot(2,1,2);imagesc(squeeze(betaBest_Lcurve(level,:,:)));title('betas');xlabel('time');ylabel('ROI')
% saveas(gcf,['figures' filesep 'betasTemplateERP_1'],'png')

% % Optimization of the regularization parameter according to the gcv error
% [u,s,v] = csvd(stackAvMap);
% %[lambda, gcvErrorMinNorm, lambdaGridMinNorm] = gcv(u,s,b,'Tikh', nLambda);
% [lambda,rho,eta,lambdaGridMinNorm] = l_curve(u,s,stackY,'Tikh');
