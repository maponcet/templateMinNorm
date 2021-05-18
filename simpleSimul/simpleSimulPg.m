clearvars;close all;

% Add folder and subfolders in path
addpath(genpath('/Volumes/Amrutam/Marlene/JUSTIN/avROImap/PlosOne/github-archive/'));
dirList = dir(['forwardAllEGI' filesep 'forward*']);
for iSubj=1:length(dirList)
    tmpName = dirList(iSubj).name(21:end);
    subjectList(iSubj) = str2num(tmpName(1:end-4));
end
numSubs = length(subjectList);

SNRlevel = 0.1; % noise level
phase = randi([2,10], 1); % for simulation. why rand 2,10 degrees... that's nothing anyway!
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
numROIs = length(listROIs);

% load average map of ROIs (128 elec x 18 ROIs)
load('averageMap50.mat') 

signalROIs = {'V1-R','MT-R','V1-L','MT-L'}; % simulated signal
% find the ROI index which corresponds to the one in signalROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(avMap.roiNames,signalROIs{x}),1:length(signalROIs),'uni',false));

%% GET ROIs and simulate signals
% do for each sbj separately but on 18 ROIs, not mesh
% so first need to get the ROI per sbj
roiMap=zeros(128,numROIs,numSubs);
roiStacked = [];
for ff=1:numSubs
    clear fwdMatrix roiInfo
    fprintf('Loading participant ROI information %d \n',ff);
    % fwd file
    load(['forwardAllEGI/forwardAndRois-skeri' num2str(subjectList(ff),'%04.f') '.mat']);
    % check that all rois are there
    if sum(ismember(listROIs, {roiInfo.name})) ~= numROIs
        fprintf('nb of ROIs mismatch for %s \n',listFiles(ff).name)
    else
        % go through each ROI and average all the mesh indexes corresponding to that ROI
        for rr=1:numROIs
            clear indexROI
            indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            [roiMap(:,rr,ff)] = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
        end
        roiStacked = blkdiag(roiStacked, roiMap(:,:,ff));
    end
end
roiStacked = bsxfun(@minus,roiStacked, mean(roiStacked));

for ss=1:numSubs
    [Y(128*(ss-1)+1:128*ss,:),source{ss},signal{ss}] = simpleSimulation(listROIs,roiMap(:,:,ss),ac_sources,SNRlevel,phase);
    % source contains the "strength" of the areas while signal is similar
    % but multiplied by y (so is over time). signal is ROIxtime without noise simulated data for each sbj
    
    % or could do with average map (= same brain for each sbj - activity is
    % different for each cluster but same cluster location)
    [Yavg(128*(ss-1)+1:128*ss,:),sourceAvg{ss},signalAvg{ss}] = simpleSimulation(listROIs,avMap.activity,ac_sources,SNRlevel,phase);
end

%  reduce dimension data (Ylo)
% =PCA denoised version of Y (denoised by truncation of the SVD)
numCols = 2; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
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

nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm 
[betaReg, ~, lambdaReg] = minimum_norm(avMap.activity, grandMeanData, nLambdaRidge);
[betaMinNorm, ~, lambdaMinNorm] = minimum_norm(roiStacked, Ylo, nLambdaRidge);
% betaMinNorm is numSubs*numROIs x time so need to extract betaMinNorm of a single sbj 
% for each ROI separately then average across sbj
regionTmp = arrayfun(@(x) betaMinNorm(1+numROIs*(x-1):numROIs*x,:),1:numSubs,'uni',false);
regionMinNorm =reshape(cell2mat(regionTmp),numROIs,size(Ylo,2),numSubs);
regionMinNormAvg = mean(regionMinNorm,3);

meanYhatReg = avMap.activity * betaReg;
YhatMN=roiStacked*betaMinNorm;
unstackedYhatMN = reshape(YhatMN,128,numSubs,size(Y,2)); 
grandMeanDataYhatMN = squeeze(mean(unstackedYhatMN,2));

% % plot average topo at a given time
% iT = 1;
% figure;
% subplot(1,3,1)
% title('simulation')
% plotContourOnScalp(grandMeanData(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off
% subplot(1,3,2)
% title('Regress')
% plotContourOnScalp(meanYhatReg(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off
% subplot(1,3,3)
% title('MinNorm')
% plotContourOnScalp(grandMeanDataYhatMN(:,iT),'skeri0044','datafiles/eegdata/')
% view(20,35)
% camproj('perspective')
% axis off

% % plot 1st sbj
% figure; plotContourOnScalp(Ylo(1:128,iT),'skeri0044','datafiles/eegdata/')
% figure; plotContourOnScalp(YhatMN(1:128,iT),'skeri0044','datafiles/eegdata/')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How can I do Group LASSO? I don't have the variables and not sure how to get them either  
% [betaVal{ii}, objValues{ii}, res] = get_solution_frobenius(X, Ytrans, betaInit, lambdaGrid(ii), alphaVal, tol, MAX_ITER, penalties, indices);




%% compute auc, mse, relative energy using average signal in rois


%%%%% AUC
tmp = zeros(numROIs,1);
tmp(ac_sources,:) = 1; 
% instead of using 1 and 0 in the matrix, use the real amount of simulated
% activity? Nooo rocArea only uses 1 and 0s
% AUC computed without normalising: it wouldn't change the results 
for nT = 1:size(Ylo,2)
    aucGpMinNormT(nT) = rocArea( abs(betaReg(:,nT)) , tmp );
    aucMinNormT(nT) = rocArea( abs(regionMinNormAvg(:,nT)) , tmp );
end
aucGpMinNorm = mean(aucGpMinNormT);
aucMinNormGp = mean(aucMinNormT); % auc calculated from the average activity (beta)

%%%%% relative energy
for nT = 1:size(Ylo,2)
    norm_betaReg(:,nT) = betaReg(:,nT) / max( abs(betaReg(:,nT)) ); % normalise estimated sources
    relEnergy(:,nT) = sum( abs( norm_betaReg(ac_sources,nT) ) ) / sum( abs(norm_betaReg(:,nT)) );
    norm_MinNorm(:,nT) = regionMinNormAvg(:,nT) / max( abs(regionMinNormAvg(:,nT)) ); % normalise estimated sources
    relEnergyMinNorm(:,nT) = sum( abs( norm_MinNorm(ac_sources,nT) ) ) / sum( abs(norm_MinNorm(:,nT)) );
end
energyGpNorm = mean(relEnergy);
energyMinNorm = mean(relEnergyMinNorm);

%%%%% mse
% average simulated signal across sbj
avSignal = mean(reshape(cell2mat(signal),numROIs,size(Ylo,2),numSubs),3);

n_avSigSource = avSignal / max(abs(avSignal)); % normalise real/simulated source
n_MSE = sum( (n_avSigSource - norm_betaReg).^2 ) / sum( (n_avSigSource).^2 ); 
mseGpNorm = mean(n_MSE);
n_MSE2 = sum( (n_avSigSource - norm_MinNorm).^2 ) / sum( (n_avSigSource).^2 ); 
mseMinNorm = mean(n_MSE2);


% plot activity per ROI over time
figure;rr=1;
set(gcf,'position',[100,100,800,1000])
for iRoi = 1:2:size(betaReg,1)
    subplot(9,2,1+2*(rr-1));hold on;
    plot(betaReg(iRoi,:)','m-','linewidth',2)
    plot(betaReg(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(betaReg)) max(max(betaReg))])
    title(listROIs(iRoi))
    subplot(9,2,2+2*(rr-1));hold on;
    plot(regionMinNormAvg(iRoi,:)','m-','linewidth',2)
    plot(regionMinNormAvg(iRoi+1,:)','b-','linewidth',2)
    ylim([min(min(regionMinNormAvg)) max(max(regionMinNormAvg))])
    rr=rr+1;
end





% % I could compute metrics per sbj (only for regular MinNorm) or on
% % average... completely different results...
% for ss=1:numSubs
%     for nT = 1:size(Ylo,2)
%         aucMinNormTSub(nT,ss) = rocArea( abs(regionMinNorm(:,nT,ss)) , tmp );
%     end
% end
% aucMinNormInd = mean(mean(aucMinNormTSub)); % auc calculated for each sbj then average
% for ss=1:numSubs
%     for nT = 1:size(Ylo,2)
%         norm_MinNormI(:,nT,ss) = regionMinNorm(:,nT,ss) / max( abs(regionMinNorm(:,nT,ss)) ); % normalise estimated sources
%         relEnergyMinNormI(nT,ss) = sum( abs( norm_MinNormI(ac_sources,nT) ) ) / sum( abs(norm_MinNormI(:,nT)) );
%     end
% end
% energyMinNormInd = mean(mean(relEnergyMinNormI));
% % not sure how to normalise - here by sbj but could have been for the
% % entire group...
% indSignal = reshape(cell2mat(signal),numROIs,size(Ylo,2),numSubs);
% for ss=1:numSubs
%     n_indSigSource = indSignal(:,:,ss) / max(abs(indSignal(:,:,ss))); % normalise real/simulated source
%     n_indMSE(ss) = mean(sum( (n_indSigSource - norm_MinNormI(:,:,ss)).^2 ) / sum( (n_indSigSource).^2 )); 
% end
% mseMinNormInd = mean(n_indMSE);
