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

allCombi = nchoosek(1:numROIs,2); % all possible combinations of areas

for simulation = 1:length(allCombi)
    fprintf('Simulation %d / %d \n',simulation,length(allCombi));
    signalROIs = listROIs(allCombi(simulation,:)); % simulated signal
    
    % find the ROI index which corresponds to the one in signalROIs
    ac_sources = cell2mat(arrayfun(@(x) cellfind(avMap.roiNames,signalROIs{x}),1:length(signalROIs),'uni',false));
    
    for repBoot =1:20
        clear Y source signal Ylo
        fprintf('bootstrap %d \n',repBoot);

        %% GET ROIs and simulate signals
        % do for each sbj separately but on 18 ROIs, not mesh
        % so first need to get the ROI per sbj
        roiMap=zeros(128,numROIs,numSubs);
        roiStacked = [];
        for ff=1:numSubs
            clear fwdMatrix roiInfo
            
            % fwd file
            load(['forwardAllEGI/forwardAndRois-skeri' num2str(subjectList(ff),'%04.f') '.mat']);
            % check that all rois are there
%             if sum(ismember(listROIs, {roiInfo.name})) ~= numROIs
%                 fprintf('nb of ROIs mismatch for %s \n',listFiles(ff).name)
%             else
                % go through each ROI and average all the mesh indexes corresponding to that ROI
                for rr=1:numROIs
                    clear indexROI
                    indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
                    [roiMap(:,rr,ff)] = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
                end
                roiStacked = blkdiag(roiStacked, roiMap(:,:,ff));
%             end
        end
        roiStacked = bsxfun(@minus,roiStacked, mean(roiStacked));
        
        for ss=1:numSubs
            [Y(128*(ss-1)+1:128*ss,:),source{ss},signal{ss}] = simpleSimulation(listROIs,roiMap(:,:,ss),ac_sources,SNRlevel,phase);
            % source contains the "strength" of the areas while signal is similar
            % but multiplied by y (so is over time). signal is ROIxtime without noise simulated data for each sbj
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
        aucGpMinNorm(repBoot,simulation) = mean(aucGpMinNormT);
        aucMinNorm(repBoot,simulation) = mean(aucMinNormT); % auc calculated from the average activity (beta)
        
        %%%%% relative energy
        for nT = 1:size(Ylo,2)
            norm_betaReg(:,nT) = betaReg(:,nT) / max( abs(betaReg(:,nT)) ); % normalise estimated sources
            relEnergy(:,nT) = sum( abs( norm_betaReg(ac_sources,nT) ) ) / sum( abs(norm_betaReg(:,nT)) );
            norm_MinNorm(:,nT) = regionMinNormAvg(:,nT) / max( abs(regionMinNormAvg(:,nT)) ); % normalise estimated sources
            relEnergyMinNorm(:,nT) = sum( abs( norm_MinNorm(ac_sources,nT) ) ) / sum( abs(norm_MinNorm(:,nT)) );
        end
        energyGpNorm(repBoot,simulation) = mean(relEnergy);
        energyMinNorm(repBoot,simulation) = mean(relEnergyMinNorm);
        
        %%%%% mse
        % average simulated signal across sbj
        avSignal = mean(reshape(cell2mat(signal),numROIs,size(Ylo,2),numSubs),3);
        
        n_avSigSource = avSignal / max(abs(avSignal)); % normalise real/simulated source
        n_MSE = sum( (n_avSigSource - norm_betaReg).^2 ) / sum( (n_avSigSource).^2 );
        mseGpNorm(repBoot,simulation) = mean(n_MSE);
        n_MSE2 = sum( (n_avSigSource - norm_MinNorm).^2 ) / sum( (n_avSigSource).^2 );
        mseMinNorm(repBoot,simulation) = mean(n_MSE2);
        
    end % end boostrap
    
end % different activated sources

% I'll have to reshape with a diagonal....  
imagesc(mean(aucGpMinNorm)); colorbar
imagesc(mean(aucMinNorm)); colorbar

