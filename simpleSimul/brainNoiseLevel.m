%%% the level of noise added to the forward varies for each participant +
%%% each simulation (ie which area(s) is(are) active)!!

% equation to find BrainNoise for a given SNR for V1-MT
% x = 1/sqrt(10*y); for step (x is noise, y is SNR) (10 comes from fitting: y= x.^-2*.1)
% x = 1/sqrt(2*y); for ERP (2 comes from fitting: y= x.^-2*.5)
% equation to find BrainNoise for V2V4: y= x.^-2*.2 ; y= x.^-2*1.2

clearvars; close all

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);


%% LOAD FWD
listSub = [44 5 11 23 38 9 42]; 
for iSub=1:length(listSub)
    pickS = listSub(iSub);
    load([dataPath dirList(pickS).name])
    fullFwd{iSub} = fwdMatrix;
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
        % save the index for each ROI
        idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
    end

    %% Simulate sources
    sourceL = {'V1-L','MT-L'};
    sourceR = {'V1-R','MT-R'};
    activeROIs = [sourceL,sourceR];
%     activeROIs = {'V1-L'};
    ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));
    % ERP & baseline timewindow
    winBase = 1:45;
    winERP = 46:90;
    % ERP source
    srcERP = zeros(numROIs,45*2); % ERP is a 45 ms window (due to timeline in the pg)
    srcERP(:,winERP) = createSourceERP(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
%     srcERP(:,winERP) = createSourceERP(numROIs,ac_sources);
    % step source
    srcStep = zeros(numROIs,45*2);
    srcStep(:,winERP) = 1;

    %%% Simulate scalp activity (Y)
    sourceData = zeros(size(fullFwd{iSub},2) , size(srcERP,2));
    sourceStep = zeros(size(fullFwd{iSub},2) , size(srcERP,2));
    for ss=1:length(ac_sources)
        sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
        sourceStep(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcStep(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
    end
    yERP = fullFwd{iSub} * sourceData;
    yStep = fullFwd{iSub} * sourceStep;


%     %%%% Add gaussian noise on ELECTRODES
%     sig = rms(rms(yStep(:,winERP))); % get amount of signal
%     SNRlevel = 0.1;
%     % check SNR using bootstrap
%     for rep = 1:1000
%         % gaussian noise around mean signal
%         noiseGauss = sig/sqrt(SNRlevel) * randn(size(yStep,1),size(yStep,2)); % should be sqrt(noise) to convert to amplitude (as is the signal)
%         yStepFull = yStep + noiseGauss;
%     
%         % ERP noise around signal variance
%         [noiseERP] = add_ERPnoise_with_SNR( yERP, SNRlevel, winERP);
%         yERPFull = yERP + noiseERP;
%     
%         stepWin = yStepFull(:,winERP); preStepWin = yStepFull(:,winBase);
%         erpWin = yERPFull(:,winERP); preWin = yERPFull(:,winBase);
%     
%         aa(rep) = (rms(yStepFull(:))/rms(noiseGauss(:)))^2 - 1; % half SNR: half win does not have signal
%         bb(rep)= (rms(stepWin(:))/rms(preStepWin(:)))^2 -1;
%         cc(rep) = (rms(yERPFull(:))/rms(noiseERP(:)))^2 - 1; % half SNR: half win does not have signal
%         dd(rep)= (rms(erpWin(:))/rms(preWin(:)))^2 -1;
%     end
%     mean(aa)
%     mean(bb)
%     mean(cc)
%     mean(dd)



    %%%% Add gaussian noise on FWD -> correlated noise
    % test multiple amount to compare with SNR level
    clear noiseGauss stepSNR erpSNR
    % sigLevel = [1:20 30 40 50 100:100:500]; % [1 10 25 50 100 250];
    noiseLevel = [0.001 0.005 0.01:0.01:0.1 0.2:0.2:1 1.5];
    for nn=1:length(noiseLevel)
        clear bb dd
        for rep = 1:100
            %     noiseGauss = 1/sigLevel(nn) * randn(size(fullFwd{iSub},2),size(yStep,2));
            noiseGauss = noiseLevel(nn) * randn(size(fullFwd{iSub},2),size(yStep,2));
            brainNoise = fullFwd{iSub} * noiseGauss; % correlated noise
            yStepBrain = yStep + brainNoise;
            yERPBrain = yERP + brainNoise;

            stepWin = yStepBrain(:,winERP); preStepWin = yStepBrain(:,winBase);
            erpWin = yERPBrain(:,winERP); preWin = yERPBrain(:,winBase);

            bb(rep)= (rms(stepWin(:))/rms(preStepWin(:)))^2 -1;
            dd(rep)= (rms(erpWin(:))/rms(preWin(:)))^2 -1;
        end
        stepSNR(nn) = mean(bb);
        erpSNR(nn) = mean(dd);
    end


    % linear change on loglog -> POWER
    % slope = Y/X
%     slope = (log(stepSNR(end)) - log(stepSNR(1))) / (log(sigLevel(end)) - log(sigLevel(1)));

    figure;scatter(noiseLevel,stepSNR)
    xlabel('brainNoise');ylabel('achieved SNR')
    set(gca,'yscale','log','xscale','log')
    hold on;
    x= linspace(0,2, 10000);
    y= x.^-2*.1; % equation to find SNR when BrainNoise is known = 1/x^2 *.1
    loglog(x,y)
    expEq = fittype('x^a*b','coeff',{'a','b'});
    [curveparam ,goodfit] = fit(noiseLevel',stepSNR',expEq,'Startpoint',[-2 .1]);
    fitSlope(iSub) = curveparam.a;
    fitInter(iSub) = curveparam.b;
    rFit(iSub) = goodfit.rsquare;

    figure;scatter(noiseLevel,erpSNR)
    xlabel('brainNoise');ylabel('achieved SNR')
    set(gca,'yscale','log','xscale','log')
    hold on;
    x= linspace(0,2, 10000);
    y= x.^-2*.5; % equation to find SNR when BrainNoise is known = 1/x^2 *.1
    loglog(x,y)
    [curveparam2,goodfit] = fit(noiseLevel',erpSNR',expEq,'Startpoint',[-2 1]);
    fitSlopeERP(iSub) = curveparam2.a;
    fitInterERP(iSub) = curveparam2.b;
    rFit(iSub) = goodfit.rsquare;

end
mean(fitSlope)
mean(fitInter)
mean(fitSlopeERP)
mean(fitInterERP)
median(fitSlope)
median(fitInter)
median(fitSlopeERP)
median(fitInterERP)


% equation to find BrainNoise for a given SNR for V1-MT
% x = 1/sqrt(10*y); for step (x is noise, y is SNR) (10 comes from fitting: y= x.^-2*.1)
% x = 1/sqrt(2*y); for ERP (2 comes from fitting: y= x.^-2*.5)


SNRlevel = 10;
clear aa bb 
for rep = 1:100

    noiseGauss = 1/sqrt(10*SNRlevel) * randn(size(fullFwd{iSub},2),size(yStep,2));
    brainNoise = fullFwd{iSub} * noiseGauss; % correlated noise
    yStepBrain = yStep + brainNoise;

    stepWin = yStepBrain(:,winERP); preStepWin = yStepBrain(:,winBase);
%     (rms(stepWin(:))/rms(preStepWin(:)))^2 -1
%     (rms(yStepBrain(:))/rms(brainNoise(:)))^2 - 1 % half SNR: half win does not have signal

    noiseGauss2 = 1/sqrt(2*SNRlevel) * randn(size(fullFwd{iSub},2),size(yERP,2));
    brainNoise2 = fullFwd{iSub} * noiseGauss2; % correlated noise
    yERPBrain = yERP + brainNoise2;

    erpWin = yERPBrain(:,winERP); preWin = yERPBrain(:,winBase);
%     (rms(erpWin(:))/rms(preWin(:)))^2 -1
%     (rms(yERPBrain(:))/rms(brainNoise2(:)))^2 - 1 % half SNR: half win does not have signal

    aa(rep)= (rms(stepWin(:))/rms(preStepWin(:)))^2 -1;
    bb(rep)= (rms(erpWin(:))/rms(preWin(:)))^2 -1;
end
mean(aa)
mean(bb)

