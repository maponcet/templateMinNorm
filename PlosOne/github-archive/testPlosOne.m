clearvars;close all;

% Y=data, G=source

% Add folder and subfolders in path
addpath(genpath('/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/'));

simulSignal=1; % look at simulation vs. real data (EEG)
if simulSignal
     % include all forwardSolutions
    dirList = dir(['datafiles' filesep 'forwardSolutions' filesep 'forward*']);
    for iSubj=1:length(dirList)
        tmpName = dirList(iSubj).name(21:end);
        subjectList(iSubj) = str2num(tmpName(1:end-4));
    end
    % missing ROIs for some sbj which creates an error when getting ROIs 
    % LOC missing for 60, 72
    subjectList = setdiff(subjectList,[60 72]);
else
   % only include EEG sbj
    dirList = dir(['datafiles' filesep 'eegdata' filesep 'skeri*']);
    for iSubj=1:length(dirList)
        subjectList(iSubj) = str2num(dirList(iSubj).name(end-3:end));
    end
end
numSubs = length(subjectList);

% main parameters
numComponents = 3; % nb of principal components that are kept?? it is for xlist and vlist...
signalROIs = {'V1-R','MT-R','V1-L','MT-L'}; % {'V2v-L','V4-R'} simulate signal coming from these ROIs (PlosOne)
SNRlevel = 0.1;
initStrct = load('Subject_48_initialization'); % not sure what this is but used in generateData
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm 
numCols = 5; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
nLambda = 30; % always end picking the max? (~120) or lasso.. 
alphaVal = 1.0817e4;
MAX_ITER = 1e6; % max iteration for group lasso solution (get_solution_frobenius)
lowPassNF1 = 0; % 0 or 1 lowPass for the EEG data (Y). Attention: will need to be consistent when I average the signal in the regression

% initialise variables
stackedForwards = [];
allSubjForwards = {};
ridgeSizes = zeros(1, numSubs);

%%%%%%%%%%%%%%%%%%
%% GET ROIs 

for s=1:numSubs
    fprintf('Loading participant ROI information %d \n',s);
    % ROI file
    structure = load(['ROI_correlation_subj_' num2str(subjectList(s))]);
    ROIs{s} = structure.ROIs;
    roiIdx{s} = ROIs{s}.ndx;
    numROIs = length(roiIdx{1}); % all subjects should have same number of ROIs
    % fwd file
    structure = load(['forwardAndRois-skeri' num2str(subjectList(s))]);
    fwdMatrix = structure.fwdMatrix;
    ridgeSizes(s) = numel(cat(2,roiIdx{s}{:})); % total ROI size by subject
    
    % stackedForwards - important but not sure why yet
    % this only keeps the ROIs whereas allSubjForwards is for whole brain. 
    stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[roiIdx{s}{:}]));
    % allSubjForwards used in the next step for generating Y
    allSubjForwards{s} = fwdMatrix;
    
    % make Xlist and Vlist
    % get the first x principle components of the part of the fwdMatrix corresponding to each ROI
    % Xlist is the list of derived forward matrices, each of dimension (n x numComponents)
    % Vlist contains the matrices used to reverse the PCA transformation
    xList{s} = cell(1,numROIs); vList{s} = cell(1,numROIs);
    [xList{s}(:),vList{s}(:)] = arrayfun(@(x) get_principal_components(fwdMatrix(:, roiIdx{s}{x}),numComponents),1:numROIs,'uni',0);
end

%%%%%%%%%%%%%%%%%%
 %% GENERATE Y 
 %% either from simulation or from the EEG data
if simulSignal
    % uses the fwdMatrix loaded for each participant in the previous step
    phase = randi([2,10], 1); % only randomize phase for first subject
%     phase = 0;
    for s = 1:numSubs
        % find the ROI index which corresponds to the one in signalROIs
        idx = cell2mat(arrayfun(@(x) cellfind(lower(ROIs{s}.name),lower(signalROIs{x})),1:length(signalROIs),'uni',false));
        if numel(idx)<numel(signalROIs) % check if ROIs exist
            error('Signal ROIs do not exist!')
        end
        % generate signal (Y is stacked: 128*numSubjects x timePoints)
        [Y(128*(s-1)+1:128*s,:), source{s}, signal{s}, noise] = GenerateData(ROIs{s},idx,initStrct.VertConn,allSubjForwards{s},SNRlevel, phase);
    end
else % Y contains stacked up data
    for s = 1:numSubs
        Axx = load(['skeri' num2str(subjectList(s),'%04.f') filesep 'Exp_MATL_HCN_128_Avg' filesep 'Axx_c008.mat']);
        % Axx_c008 is the condition used for plotting in PlosOne paper
        if lowPassNF1
            harm2keep = 10;
            nf1 = Axx.i1F1;
            axxIdx = (nf1:nf1:harm2keep*nf1)+1;
            dftIdx = (1:harm2keep)+1;
            dft = dftmtx(size(Axx.Wave,1));
            sinWaveforms = imag(dft);
            cosWaveforms = real(dft);
            wave = cosWaveforms(:,dftIdx)*Axx.Cos(axxIdx,:)-sinWaveforms(:,dftIdx)*Axx.Sin(axxIdx,:);
            Axx.Wave = wave;
        end
        % Y is 128*9 sbj x 780 timepoints
        Y(size(Axx.Wave,2)*(s-1)+1:size(Axx.Wave,2)*s,:) = Axx.Wave';
    end
end

%  reduce dimension data (Ylo)
% =PCA denoised version of Y (denoised by truncation of the SVD)
  


%%%%%%%%%%%%%%%%%%
%%  GENERATE X AND V
disp('Generate X and V');
X = []; V = [];
for g = 1:numROIs
    tempX = [];
    tempV = [];
    for s = 1:numSubs
        tempX = blkdiag(tempX, xList{s}{g});
        tempV = blkdiag(tempV, vList{s}{g}(:,1:numComponents));
    end
    X = [X, tempX];
    V = blkdiag(V, tempV);
end
grpSizes = numComponents*numSubs*ones(1,numROIs);
indices = get_indices(grpSizes);
penalties = get_group_penalties(X, indices);


%%%%%%%%%%%%%%%%%%
% done with generating Y and X
% center Y, X, stackedForwards
Y2 = bsxfun(@minus,Y, mean(Y));
X = bsxfun(@minus,X, mean(X));
stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));




%% Create data topographies from each participant.
unstackedData = reshape(Ylo,128,numSubs,size(Ylo,2)); % 128electrodes x 9subs x 91 timepoints
grandMeanData = squeeze(mean(unstackedData,2));
if simulSignal
    iT = 1;
else
    iT = 196; %Time index to plot. 196th sample = 251.3 ms (PlosOne)
end

if simulSignal
    figure;
    plotContourOnScalp(grandMeanData(:,iT),'skeri0044','datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
else
    for iSubj = 1:numSubs
        figure;
        thisSubj = dirList(iSubj).name;
        %         subplot(1,9,iSubj)
        plotContourOnScalp(unstackedData(:,iSubj,iT),thisSubj,'datafiles/eegdata/')
        view(20,35)
        camproj('perspective')
        axis off
    end
end



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Minimum Norm Solution - takes a LONG time


if simulSignal
    % % traditional minimum norm solution
    % [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Y, nLambdaRidge);
    % lambdaMinNorm = lambdaMinNorm^2;
    % rsquaredMinNorm = 1 - (norm(Y-stackedForwards*betaMinNorm, 'fro')^2/n) / ssTotal;
    
    [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Ylo, nLambdaRidge);
    lambdaMinNorm = lambdaMinNorm^2;
    rsquaredMinNorm = 1 - (norm(Y-stackedForwards*betaMinNorm, 'fro')^2/n) / ssTotalLo;
    
    save('simul_calculatedMinNorm.mat','betaMinNorm','lambdaMinNorm','rsquaredMinNorm');
else
    load('precalculatedMinNorm.mat')
end
YhatMN=stackedForwards*betaMinNorm;
unstackedYhatMN = reshape(YhatMN,128,numSubs,size(Y,2)); % 128electrodes x 9subs x 91 timepoints
grandMeanDataYhatMN = squeeze(mean(unstackedYhatMN,2));
figure;
title('min norm')
plotContourOnScalp(grandMeanDataYhatMN(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Group LASSO
% use first 2 columns of v as time basis
% PK: moved this down to avoid conflict with v generated above

disp('Starting group LASSO');

disp('Reducing dimensionality of data');
%use first 2/5 columns of v as time basis
[~, ~, v] = svd(Y);
Ytrans = Y * v(:, 1:numCols);
% center
Ytrans = bsxfun(@minus,Ytrans, mean(Ytrans));

% sequence of lambda values
lambdaMax = max(cell2mat(arrayfun(@(x) norm(X(:,indices{x})'*Ytrans, 'fro')/penalties(x),1:numROIs,'uni',0)));
lambdaMax = lambdaMax + 1e-4;
lambdaGrid = lambdaMax * (0.01.^(0:1/(nLambda-1):1));
tol = min(penalties) * lambdaGrid(end) * 1e-5;
if alphaVal > 0
    tol = min([tol, 2*alphaVal*1e-5]);
end

ridgeRange = [0 cumsum(ridgeSizes)];
roiSizes = zeros(1,numROIs); %total size of each region summed over all subjects
for g = 1:numROIs
    roiSizes(g) = sum(cell2mat(arrayfun(@(x) numel(roiIdx{x}{g}),1:numSubs,'uni',0)));
end

% OLS FIT optimal least square
betaOls = (X'*X + alphaVal*eye(size(X,2))) \ (X'*Ytrans);

% FITTING
disp(['Findig Group-LASSO solution for ' num2str(nLambda) ' regularization (lambda) values']) ;

%Find the optimal regularization parameter for the group-lasso solutions
betaInit = zeros(size(X,2), numCols);
betaVal = cell(1, nLambda);
objValues = cell(1, nLambda);
gcvError = zeros(1, nLambda);
df = zeros(1, nLambda);
indexer = cell2mat(arrayfun(@(x) return_index(roiSizes, roiIdx, x), 1:numSubs,'uni',0));
for ii = 1:nLambda
    fprintf('lambda = %d ',ii);
    [betaVal{ii}, objValues{ii}, res] = get_solution_frobenius(X, Ytrans, betaInit, lambdaGrid(ii), alphaVal, tol, MAX_ITER, penalties, indices);
    betaInit = betaVal{ii};
    betaVal{ii} = V * betaVal{ii} * v(:,1:numCols)'; %transform back to original space (permuted forward matrices)
    rss = norm(Y-stackedForwards*betaVal{ii}(indexer, :), 'fro')^2 / n;
    [gcvError(ii), df(ii)] = compute_gcv(rss, betaInit, betaOls, grpSizes, n);
end
[~, bestIndex] = min(gcvError);
YhatLASSO = stackedForwards*betaVal{bestIndex}(indexer, :);

unstackedYhatLASSO = reshape(YhatLASSO,128,numSubs,size(Y,2)); % 128electrodes x 9subs x 91 timepoints
grandMeanDataYhatLASSO = squeeze(mean(unstackedYhatLASSO,2));
figure;
title('lasso')
plotContourOnScalp(grandMeanDataYhatLASSO(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regression: Data = ROI * Betas. Matlab: Y = X*B, B = regress(Y,X)
% could use mldivide as well or B= X\Y
% get the average roi maps and regress with grandMeanData
% to add new skeri in the model can run makeManyForwards.m?
load('averageMap.mat') % 128 electrodes 18 ROIs

% % not sure yet about normalisation. Would need if different amount of noise
% % per sbj?
% % do a normalisation for each ROI -> unit norming
% % so that the betas represent microVolts (instead of microVolts/area size
% % as it is now)
% % unit norming is: all electrodes are squared and summed. These values are
% % then divided so that the total of the electrodes for each ROI (power) is
% % equal to 1
% normTerm = sqrt(sum(avMap.activity.^2,1));
% normModel = bsxfun(@rdivide,avMap.activity,normTerm);

%%%% previous code with regress. Instead use min-norm function which uses a
%%%% regularisation parameter
% betaRegress=zeros(length(avMap.roiNames),size(grandMeanData,2));
% meanYhatRegress=zeros(size(grandMeanData));
% for nT=1:size(grandMeanData,2)
%     betaRegress(:,nT) = regress(grandMeanData(:,nT), avMap.activity);
%     meanYhatRegress(:,nT) = avMap.activity * betaRegress(:,nT);
% end
% % indexR = 2:2:18;
% % for nT=1:size(grandMeanData,2)
% %     betaRegress(:,nT) = regress(grandMeanData(:,nT), avMap.activity(:,indexR));
% %     meanYhatRegress(:,nT) = avMap.activity(:,indexR) * betaRegress(:,nT);
% % end
% % for nT=1:size(grandMeanData,2)
% %     betaRegressNorm(:,nT) = regress(grandMeanData(:,nT), normModel(:,indexR));
% %     meanYhatRegressNorm(:,nT) = normModel(:,indexR) * betaRegress(:,nT);
% % end
% iT=25;
% figure;
% title('regression')
% plotContourOnScalp(meanYhatRegress(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off


figure;
for iRoi = 1:size(betaRegress,1)
    subplot(9,2,iRoi);hold on;
    plot(betaRegress(iRoi,:)','k-','linewidth',3)
    ylim([min(min(betaRegress)) max(max(betaRegress))])
    title(avMap.roiNames(iRoi))
end
%   figure;
% for iRoi = 1:size(betaRegress,1)
%     subplot(9,1,iRoi);hold on;
%     plot(betaRegressNorm(iRoi,:)','k-','linewidth',3)
%     ylim([min(min(betaRegressNorm)) max(max(betaRegressNorm))])
%     title(avMap.roiNames(indexR(iRoi)))
% end



[betaReg2, ~, lambdaReg] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
meanYhatReg2 = avMap.activity * betaReg2;
figure;
title('regression')
plotContourOnScalp(meanYhatReg2(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off
figure;
for iRoi = 1:size(betaReg2,1)
    subplot(9,2,iRoi);hold on;
    plot(betaReg2(iRoi,:)','k-','linewidth',3)
    ylim([min(min(betaReg2)) max(max(betaReg2))])
    title(avMap.roiNames(iRoi))
end
set(gcf,'position',[100,100,500,1000])

% No idea if the following lines are correct
% nMean = numel(grandMeanData);
% ssMeanLo = norm(grandMeanData, 'fro')^2 / nMean;
% lambdaReg = lambdaReg^2;
% rsquaredReg = 1 - (norm(grandMeanData-meanYhatReg2, 'fro')^2/n) / ssMeanLo;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
% all topo on the same graph
iT = 25;
figure;
subplot(2,2,1)
title('simulation T=25')
plotContourOnScalp(grandMeanData(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off
subplot(2,2,2)
title('min norm')
plotContourOnScalp(grandMeanDataYhatMN(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off
subplot(2,2,3)
title('lasso')
plotContourOnScalp(grandMeanDataYhatLASSO(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off
subplot(2,2,4)
title('regression')
plotContourOnScalp(meanYhatRegress(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off








%% compute average region activity
%%% average activity per ROI across subj for minNorm and lasso solutions
regionActivityMinNorm = zeros(numROIs, size(Y, 2),numSubs);
regionActivityLasso = zeros(numROIs, size(Y, 2),numSubs);       
for s = 1:numSubs
    % correct indexing of rois for each subject s
    range = cell2mat(arrayfun(@(x)  numel(roiIdx{s}{x}),1:numROIs,'uni',false));
    % does the same as line above
%     for i = 1:numROIs
%         range(i) = numel(roiIdx{s}{i}); 
%     end
    range = [0 cumsum(range)];
    
    tempMinNorm = betaMinNorm(ridgeRange(s)+1:ridgeRange(s+1), :);
    temp = betaVal{bestIndex}(return_index(roiSizes, roiIdx, s), :);
    regionActivityMinNorm(:,:,s) = cell2mat(arrayfun(@(x) mean(tempMinNorm(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
    regionActivityLasso(:,:,s) = cell2mat(arrayfun(@(x) mean(temp(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
end
meanRegionActivityMinNorm = mean(regionActivityMinNorm,3);
meanRegionActivityLasso = mean(regionActivityLasso,3);

% this is for real data with a baseline???
% baselineIdx = 700:780; 
% meanRegionActivityLasso2 = bsxfun(@minus,meanRegionActivityLasso,mean(meanRegionActivityLasso(:,baselineIdx),2));
% meanRegionActivityMinNorm2 = bsxfun(@minus,regionActivityMinNorm,mean(meanRegionActivityMinNorm(:,baselineIdx),2));

%% plot ROI activity over 'time'
rightIdx = cell2mat(arrayfun(@(x) ~isempty(strfind(ROIs{1}.name{x},'-R')), 1:length(ROIs{1}.name),'uni',false));
selectedRois = find(rightIdx==0);% select right==1 or left==0 hemisphere rois
nRoi=length(selectedRois);
figure;
for iRoi = 1:nRoi
    subplot(nRoi,2,iRoi+(iRoi-1));hold on;
    plot(meanRegionActivityMinNorm(selectedRois(iRoi),:)','k-','linewidth',3)
%     ylim([-2 2])
    ylim([min(min(meanRegionActivityMinNorm)) max(max(meanRegionActivityMinNorm))])
    title(ROIs{1}.name{selectedRois(iRoi)})
    subplot(nRoi,2,iRoi+iRoi);hold on;
    plot(meanRegionActivityLasso(selectedRois(iRoi),:)','k-','linewidth',3)
    ylim([min(min(meanRegionActivityLasso)) max(max(meanRegionActivityLasso))])
    title(ROIs{1}.name{selectedRois(iRoi)})
end

figure;
for iRoi = 1:size(betaRegress,1)
    subplot(9,2,iRoi);hold on;
    plot(betaRegress(iRoi,:)','k-','linewidth',3)
    ylim([min(min(betaRegress)) max(max(betaRegress))])
    title([avMap.roiNames{iRoi} ' reg'])
end
set(gcf,'position',[100,100,500,1000])
figure;
for iRoi = 1:size(meanRegionActivityMinNorm,1)
    subplot(9,2,iRoi);hold on;
    plot(meanRegionActivityMinNorm(iRoi,:)','k-','linewidth',3)
    ylim([min(min(meanRegionActivityMinNorm)) max(max(meanRegionActivityMinNorm))])
    title([ROIs{1}.name{iRoi} ' minNorm'])
end
set(gcf,'position',[100,100,500,1000])
figure;
for iRoi = 1:size(meanRegionActivityLasso,1)
    subplot(9,2,iRoi);hold on;
    plot(meanRegionActivityLasso(iRoi,:)','k-','linewidth',3)
    ylim([min(min(meanRegionActivityLasso)) max(max(meanRegionActivityLasso))])
    title([ROIs{1}.name{iRoi} ' lasso'])
end
set(gcf,'position',[100,100,500,1000])


%%% Figure showing time courses from minimum-norm solution for every ROI on single plot.
figure;
subplot(2,2,1);
plot(meanRegionActivityMinNorm(rightIdx==0,:)','linewidth',2)
legend(ROIs{1}.name(rightIdx==0),'location','southeastoutside');
title('Min Norm Solution - Left Hemisphere');
ylabel('Current Source Density')
subplot(2,2,2);hold on;
plot(meanRegionActivityMinNorm(rightIdx,:)','linewidth',2)
legend(ROIs{1}.name(rightIdx),'location','southeastoutside');
title('Min Norm Solution - Right Hemisphere');
ylabel('Current Source Density')
subplot(2,2,3);
plot(meanRegionActivityLasso(rightIdx==0,:)','linewidth',2)
legend(ROIs{1}.name(rightIdx==0),'location','southeastoutside');
title('Gp lasso - Left Hemisphere');
ylabel('Current Source Density')
subplot(2,2,4);hold on;
plot(meanRegionActivityLasso(rightIdx,:)','linewidth',2)
legend(ROIs{1}.name(rightIdx),'location','southeastoutside');
title('Gp lasso - Right Hemisphere');
ylabel('Current Source Density')









%% compute auc, mse energy for min norm and lasso
% AUC quantifies how well the estimated currents detect true sources and reject false positives
% An AUC close to 1 means that the model separates the active and nonactive sets of sources well. 
% However, in our simulations, the number of inactive sources is very much larger than the number of active sources.
% estimates close and far AUC so In the end, a global AUC value can be computed as an average of these two scalars:
% AUC=1/2(AUCclose + AUCfar)
aucCloseMN = zeros(1, numSubs);
aucFarMN = zeros(1, numSubs);
mseMN = zeros(1, numSubs);
energyMN = zeros(1, numSubs);
aucCloseLASSO = zeros(1, numSubs);
aucFarLASSO = zeros(1, numSubs);
mseLASSO = zeros(1, numSubs);
energyLASSO = zeros(1, numSubs);
for ss = 1:numSubs
    [~, aucCloseMN(ss), aucFarMN(ss), mseMN(ss), ~, energyMN(ss)] = get_metrics(betaMinNorm(ridgeRange(ss)+1:ridgeRange(ss+1), :), signal{ss}, initStrct.VertConn, roiIdx{ss});
    aucGlobalMN(ss) = 1/2 * (aucCloseMN(ss) + aucFarMN(ss));
    [~, aucCloseLASSO(ss), aucFarLASSO(ss), mseLASSO(ss), ~, energyLASSO(ss)] = get_metrics(betaVal{bestIndex}(return_index(roiSizes, roiIdx, ss), :), signal{ss}, initStrct.VertConn, roiIdx{ss});    
    aucGlobalLASSO(ss) = 1/2 * (aucCloseLASSO(ss) + aucFarLASSO(ss));
end
