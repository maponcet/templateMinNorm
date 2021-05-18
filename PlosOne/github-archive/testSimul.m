% simulate activity from the average ROI (->Y)
% reduce dimension data (Ylo)
% find ind min_norm, groupLasso, group min_norm

clearvars;close all;

% Add folder and subfolders in path
addpath(genpath('/Volumes/Amrutam/Marlene/JUSTIN/avROImap/PlosOne/github-archive/'));

dirList = dir(['forwardAllEGI' filesep 'forward*']);

for iSubj=1:length(dirList)
    tmpName = dirList(iSubj).name(21:end);
    subjectList(iSubj) = str2num(tmpName(1:end-4));
end
numSubs = length(subjectList);

% main parameters
numComponents = 3; % nb of principal components that are kept?? it is for xlist and vlist...
signalROIs = {'V1-R','MT-R','V1-L','MT-L'}; % {'V2v-L','V4-R'}; % simulate signal coming from these ROIs
SNRlevel = 0.1;
initStrct = load('Subject_48_initialization'); % not sure what this is but used in generateData: VertConn=identity matrix size nbROI x nbROI
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm 
numCols = 5; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
nLambda = 30; % falways end picking the max? (~120) or lasso.. 
alphaVal = 1.0817e4;
MAX_ITER = 1e6; % max ieration for group lasso solution (get_solution_frobenius)

% initialise variables
stackedForwards = [];
allSubjForwards = {};
ridgeSizes = zeros(1, numSubs);

%%%%%%%%%%%%%%%%%%
%% GET ROIs 
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
numROIs = length(listROIs);
% go through each ROI and average all the mesh indexes corresponding to that ROI
    
for s=1:numSubs
    fprintf('Loading participant ROI information %d \n',s);
    % fwd file
    structure = load(['forwardAllEGI/forwardAndRois-skeri' num2str(subjectList(s),'%04.f') '.mat']);
    % go through each ROI to find its corresponding indexes in fwdsolution
    for rr=1:length(listROIs)
        indexROI(rr) = find(strcmp(listROIs(rr),{structure.roiInfo.name}));
    end
    ROIs{s}.ndx = {structure.roiInfo(indexROI).meshIndices};
    ROIs{s}.name = listROIs;
    roiIdx{s} = ROIs{s}.ndx;
    
    fwdMatrix = structure.fwdMatrix;
    ridgeSizes(s) = numel(cat(2,roiIdx{s}{:})); % total ROI size by subject
    
    % stackedForwards 
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
% uses the fwdMatrix loaded for each participant in the previous step
phase = randi([2,10], 1); % only randomize phase for first subject
%     phase = 0;
for s = 1:numSubs
    % find the ROI index which corresponds to the one in signalROIs
    ac_sources = cell2mat(arrayfun(@(x) cellfind(ROIs{s}.name,signalROIs{x}),1:length(signalROIs),'uni',false));
    if numel(ac_sources)<numel(signalROIs) % check if ROIs exist
        error('Signal ROIs do not exist!')
    end
    % generate signal (Y is stacked: 128*numSubjects x timePoints)
    % if very far point then not activated even with cluster_size=100
    [Y(128*(s-1)+1:128*s,:), source{s}, signal{s}, noise] = GenerateData(ROIs{s},ac_sources,initStrct.VertConn,allSubjForwards{s},SNRlevel, phase);
end

%  reduce dimension data (Ylo)
% =PCA denoised version of Y (denoised by truncation of the SVD)
n = numel(Y);
[u1, s1, v1] = svd(Y);
ssTotal = norm(Y, 'fro')^2 / n;
Ylo = u1(:,1:numCols)*s1(1:numCols,1:numCols)*v1(:, 1:numCols)';
ssTotalLo = norm(Ylo, 'fro')^2 / n;
        

unstackedData = reshape(Ylo,128,numSubs,size(Ylo,2)); 
grandMeanData = squeeze(mean(unstackedData,2));
%%%% plot Y
% I don't know what is the timecourse - 1 oscillation
figure;plot(grandMeanData','k')
ylabel('Potential (microvolts)')
% plot average 
iT = 25; %Time index to plot 1=max sinewave when phase = 1?
figure;
title('simulation')
plotContourOnScalp(grandMeanData(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off



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
Y = bsxfun(@minus,Y, mean(Y));
X = bsxfun(@minus,X, mean(X));
stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Find Minimum Norm Solution - takes a LONG time
% % traditional minimum norm solution
% [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Y, nLambdaRidge);
% lambdaMinNorm = lambdaMinNorm^2;
% rsquaredMinNorm = 1 - (norm(Y-stackedForwards*betaMinNorm, 'fro')^2/n) / ssTotal;

[betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Ylo, nLambdaRidge);
lambdaMinNorm = lambdaMinNorm^2;
rsquaredMinNorm = 1 - (norm(Y-stackedForwards*betaMinNorm, 'fro')^2/n) / ssTotalLo;
     
YhatMN=stackedForwards*betaMinNorm;
% save('simul_calculatedMinNorm50-V1MT.mat','betaMinNorm','lambdaMinNorm','rsquaredMinNorm','YhatMN');

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
save('simul_calculatedLasso50-V1MT.mat','bestIndex','betaVal','gcvError','YhatLASSO');

unstackedYhatLASSO = reshape(YhatLASSO,128,numSubs,size(Y,2)); % 128electrodes x 9subs x 91 timepoints
grandMeanDataYhatLASSO = squeeze(mean(unstackedYhatLASSO,2));
figure;
title('lasso')
plotContourOnScalp(grandMeanDataYhatLASSO(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off


%% compute average activity per region for minNorm and lasso
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
    
    tempMinNorm = betaMinNorm(ridgeRange(s)+1:ridgeRange(s+1), :); % extract sub
    temp = betaVal{bestIndex}(return_index(roiSizes, roiIdx, s), :);
    regionActivityMinNorm(:,:,s) = cell2mat(arrayfun(@(x) mean(tempMinNorm(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)'); % compute per roi
    regionActivityLasso(:,:,s) = cell2mat(arrayfun(@(x) mean(temp(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
end
meanRegionActivityMinNorm = mean(regionActivityMinNorm,3);
meanRegionActivityLasso = mean(regionActivityLasso,3);

% this is for real data with a baseline???
% baselineIdx = 700:780; 
% meanRegionActivityLasso2 = bsxfun(@minus,meanRegionActivityLasso,mean(meanRegionActivityLasso(:,baselineIdx),2));
% meanRegionActivityMinNorm2 = bsxfun(@minus,regionActivityMinNorm,mean(meanRegionActivityMinNorm(:,baselineIdx),2));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% regression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regression: Data = ROI * Betas. Matlab: Y = X*B, B = regress(Y,X)
% could use mldivide as well or B= X\Y
% Use Ylo to fit the data
% NO! Use minimum_norm function instead of regress so that a regularisation 
% parameter is included as it is the case for the other methods
load('averageMap50.mat') % 128 elec x 18 ROIs

%%%%TEST
[betaO, betaMinNormO, lambdaO,gcvErrorMinNormO, lambdaMinNorm] = minimum_norm(stackedForwards, Ylo, nLambdaRidge);

[betaReg1, betaMinNorm1, lambda1, gcvErrorMinNorm1, lambdaReg1] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
[betaReg2, betaMinNorm2, lambda2, gcvErrorMinNorm2, lambdaReg2] = minimum_norm(repmat(avMap.activity,10,1), Ylo, nLambdaRidge);
[betaReg3, betaMinNorm3, lambda3, gcvErrorMinNorm3, lambdaReg3] = minimum_norm(avMap.activity, grandMeanData, 1000);
[betaReg4, betaMinNorm4, lambda4, gcvErrorMinNorm4, lambdaReg4] = minimum_norm(stackedForwards,repmat(grandMeanData,10,1), nLambdaRidge);
figure;rr=1;
for iRoi = 1:2:size(betaReg,1)
    subplot(9,3,1+3*(rr-1));hold on;
    plot(betaReg(iRoi,:)','m-','linewidth',2)
    plot(betaReg(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(betaReg)) max(max(betaReg))])
    title(avMap.roiNames(iRoi))
    subplot(9,3,2+3*(rr-1));hold on;
    plot(betaReg3(iRoi,:)','m-','linewidth',2)
    plot(betaReg3(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(betaReg3)) max(max(betaReg3))])
    subplot(9,3,3+3*(rr-1));hold on;
    plot(betaReg4(iRoi,:)','m-','linewidth',2)
    plot(betaReg4(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(betaReg4)) max(max(betaReg4))])
    rr=rr+1;
end

%%%%ENDTEST



[betaReg, ~, lambdaReg] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
meanYhatReg = avMap.activity * betaReg;
% betaReg corresponds to mean region activity over time?? It's weights, not
% activity per se... 

figure;
title('regression')
plotContourOnScalp(meanYhatReg(:,iT),'skeri0044','datafiles/eegdata/')
view(20,35)
camproj('perspective')
axis off


% % all topo on the same graph
% iT = 1;
% figure;
% subplot(2,2,1)
% title('simulation')
% plotContourOnScalp(grandMeanData(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off
% subplot(2,2,2)
% title('min norm')
% plotContourOnScalp(grandMeanDataYhatMN(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off
% subplot(2,2,3)
% title('lasso')
% plotContourOnScalp(grandMeanDataYhatLASSO(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off
% subplot(2,2,4)
% title('regression')
% plotContourOnScalp(meanYhatReg(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off




%%% 
figure;rr=1;
for iRoi = 1:2:size(betaReg,1)
    subplot(9,3,1+3*(rr-1));hold on;
    plot(betaReg(iRoi,:)','m-','linewidth',2)
    plot(betaReg(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(betaReg)) max(max(betaReg))])
    title(avMap.roiNames(iRoi))
    subplot(9,3,2+3*(rr-1));hold on;
    plot(meanRegionActivityMinNorm(iRoi,:)','m-','linewidth',2)
    plot(meanRegionActivityMinNorm(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(meanRegionActivityMinNorm)) max(max(meanRegionActivityMinNorm))])
    subplot(9,3,3+3*(rr-1));hold on;
    plot(meanRegionActivityLasso(iRoi,:)','m-','linewidth',2)
    plot(meanRegionActivityLasso(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(meanRegionActivityLasso)) max(max(meanRegionActivityLasso))])
    rr=rr+1;
end
set(gcf,'position',[100,100,800,1000])
saveas(gcf,'figures/simulSolV1MT.png','png')


% % do a normalisation for each ROI -> unit norming
% % so that the betas represent microVolts (instead of microVolts/area size
% % as it is now)
% % unit norming is: all electrodes are squared and summed. These values are
% % then divided so that the total of the electrodes for each ROI (power) is
% % equal to 1
% modelNorm = sqrt(sum(avMap.activity.^2,1));
% model = bsxfun(@rdivide,avMap.activity,modelNorm);
% figure; imagesc(avMap.activity)
% figure; imagesc(model)
% [betaReg2, ~, lambdaReg2] = minimum_norm(model, grandMeanData, nLambdaRidge);
% meanYhatReg2 = model * betaReg2;





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


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% compute the average real (simulated) source distribution
for ss = 1:numSubs
    for iROI=1:length(listROIs)
        idx=roiIdx{ss}{iROI};
        sigSource(iROI,ss,:) = mean(signal{ss}(idx,:));
        simSource(iROI,ss) = mean(source{ss}(idx));
    end
end
avSigSource = squeeze(mean(sigSource,2));
figure;rr=1;
for iRoi = 1:2:size(betaReg,1)
    subplot(9,1,rr);hold on;
    plot(avSigSource(iRoi,:)','m-','linewidth',2)
    plot(avSigSource(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(avSigSource)) max(max(avSigSource))])
    title(listROIs(iRoi))
    rr=rr+1;
end

nameBar = categorical(listROIs);
nameBar = reordercats(nameBar,listROIs);
figure;bar(nameBar,mean(simSource,2));
ylabel('source activation')
title('simulated sources')

%% compute auc, mse, relative energy using average signal in rois
%%% 
% find the ROI index which corresponds to the one in signalROIs
% ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,signalROIs{x}),1:length(signalROIs),'uni',false));  % ndx to the activated sources
% rocArea works only for a vector
tmp = zeros(1,size(avSigSource,1));
tmp(ac_sources,:) = 1; 
% instead of using 1 and 0 in the matrix, use the real amount of simulated
% activity?
% AUC computed without normalising: it wouldn't change the results 
for nT = 1:length(avSigSource)
    aucGpMinNormT(nT) = rocArea( abs(betaReg(:,nT)) , tmp );
    aucMinNormT(nT) = rocArea( abs(meanRegionActivityMinNorm(:,nT)) , tmp );
    aucLassoT(nT) = rocArea( abs(meanRegionActivityLasso(:,nT)) , tmp );
end
aucGpMinNorm = mean(aucGpMinNormT);
aucMinNorm = mean(aucMinNormT);
aucLasso = mean(aucLassoT);

for nT = 1:length(avSigSource)
    norm_betaReg(:,nT) = betaReg(:,nT) / max( abs(betaReg(:,nT)) ); % normalise estimated sources
    relEnergy(:,nT) = sum( abs( norm_betaReg(ac_sources,nT) ) ) / sum( abs(norm_betaReg(:,nT)) );
    norm_MinNorm(:,nT) = meanRegionActivityMinNorm(:,nT) / max( abs(meanRegionActivityMinNorm(:,nT)) ); % normalise estimated sources
    relEnergyMinNorm(:,nT) = sum( abs( norm_MinNorm(ac_sources,nT) ) ) / sum( abs(norm_MinNorm(:,nT)) );
    norm_lasso(:,nT) = meanRegionActivityLasso(:,nT) / max( abs(meanRegionActivityLasso(:,nT)) ); % normalise estimated sources
    relEnergyLasso(:,nT) = sum( abs( norm_lasso(ac_sources,nT) ) ) / sum( abs(norm_lasso(:,nT)) );
end
energyGpNorm = mean(relEnergy);
energyMinNorm = mean(relEnergyMinNorm);
energyLasso = mean(relEnergyLasso);


mse = norm(avSigSource-betaReg, 'fro')^2 / numel(betaReg); %%% give WEIRD nb!

n_avSigSource = avSigSource / max(abs(avSigSource)); % normalise real/simulated source
n_MSE = sum( (n_avSigSource - norm_betaReg).^2 ) / sum( (n_avSigSource).^2 ); 
mseGpNorm = mean(n_MSE);
n_MSE2 = sum( (n_avSigSource - norm_MinNorm).^2 ) / sum( (n_avSigSource).^2 ); 
mseMinNorm = mean(n_MSE2);
n_MSE3 = sum( (n_avSigSource - norm_lasso).^2 ) / sum( (n_avSigSource).^2 ); 
mseLasso = mean(n_MSE3);

figure;bar([mean(aucGpMinNorm) mean(aucMinNorm) mean(aucLasso)])
title('AUC')
figure;bar([mean(mseGpNorm) mean(mseMinNorm) mean(mseLasso)])
title('MSE')
figure;bar([mean(energyGpNorm) mean(energyMinNorm) mean(energyLasso)])
title('Energy')


% %%%%% TEST
% [bY, bSource, bSignal, noise] = GenerateData(ROIs{1},1,initStrct.VertConn,allSubjForwards{1},SNRlevel, phase);
% bSource(roiIdx{1}{1})
% bSource(roiIdx{1}{2})
% bSource(roiIdx{1}{3})
% 
% for s=1:2
% [testY, testSource{s}, testSignal{s}, noise] = GenerateData(ROIs{s},[2 18 1 17],initStrct.VertConn,allSubjForwards{s},SNRlevel, phase);
% end
% testSource{s}(roiIdx{s}{18})
% testSignal(roiIdx{s}{18},:)
% 
% mean(testSource{s}(roiIdx{s}{1}))
% mean(testSource{s}(roiIdx{s}{3}))

