clearvars;close all;
% Retrieve sources using real data from Lim 2017 using the Lasso method

lowPassNF1 = 1; % filter 1 or not 0
numFq2keep = 5; % nb of harmonics to keep in the signal
totBoot = 500; %500 or 0 for no bootstap
numComponentsEEG = 5; % Nb of principal components for the EEG
numComponents = 3; % Nb of principal components for the group lasso
addpath(genpath([pwd filesep 'subfunctions']))
nLambdaRidge = 50;
nLambda = 30;
alphaVal = 1.0817e4;
MAX_ITER = 1e6;

sbjList = dir('realData/eegdata/skeri*');
numSubs = length(sbjList);
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
numROIs = length(listROIs);

%% start bootstrap
for bb=1:totBoot
    
    if mod(bb,10)==0;fprintf('boot nb %d',bb);end
    pickN = randi(numSubs,1,numSubs);
    stackedForwards = [];fwdMatrix = [];fullFwd={};roiIdx={};xList={};vList={};

%% Load forwards
for iSub=1:numSubs
    % fwd file
    load(['realData/forwardSolutions/forwardAndRois-' sbjList(pickN(iSub)).name '.mat'])
    fullFwd{iSub} = fwdMatrix;
    % go through each ROI and save the corresponding fwdMesh values
    % corresponding to the indexes of that ROI
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiIdx{iSub}{:,rr} = roiInfo(indexROI).meshIndices;
    end
    ridgeSizes(iSub) = numel(cat(2,roiIdx{iSub}{:})); % total ROI size by subject
    stackedForwards = blkdiag(stackedForwards, fwdMatrix(:,[roiIdx{iSub}{:}]));
    % make Xlist and Vlist
    % get the first x principle components of the part of the fwdMatrix corresponding to each ROI
    xList{iSub} = cell(1,numROIs); vList{iSub} = cell(1,numROIs);
    [xList{iSub}(:),vList{iSub}(:)] = arrayfun(@(x) get_principal_components(fwdMatrix(:, roiIdx{iSub}{x}),numComponents),1:numROIs,'uni',0);
end
    
 %% Load EEG data    
clear Y;
for iSub=1:numSubs
    % eegFile, Lim uses condNmbr = 8;
    Axx = load(['realData/eegdata/' sbjList(pickN(iSub)).name '/Exp_MATL_HCN_128_Avg/Axx_c008.mat']);
    
    if lowPassNF1
        nf1 = Axx.i1F1;
        axxIdx = (nf1:nf1:numFq2keep*nf1)+1; 
        dftIdx = (1:numFq2keep)+1; % dftIdx has 1 step per Hz whereas axxIdx is 0.5Hz 
        dft = dftmtx(size(Axx.Wave,1)); 
        sinWaveforms = imag(dft);
        cosWaveforms = real(dft);
        wave = cosWaveforms(:,dftIdx)*Axx.Cos(axxIdx,:)-sinWaveforms(:,dftIdx)*Axx.Sin(axxIdx,:);
        Axx.Wave = wave;  
    end
    Y(iSub,:,:) = Axx.Wave';
end
t = 0:Axx.dTms:(Axx.nT-1)*Axx.dTms;
        

%%  GENERATE X AND V
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

%% "PCA"
tempY = permute(Y,[2 1 3]);
stackY = reshape(tempY,[size(Y,2)*numSubs,size(Y,3)]);
% center Y, X, stackedForwards
stackY = bsxfun(@minus,stackY, mean(stackY));
X = bsxfun(@minus,X, mean(X));
stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));
% PCA
[u1, s1, v1] = svd(stackY);
Ylo = u1(:,1:numComponentsEEG)*s1(1:numComponentsEEG,1:numComponentsEEG)*v1(:, 1:numComponentsEEG)';
n = numel(stackY);
ssTotal = norm(Ylo, 'fro')^2 / n;

%% Lasso parameters
% disp('Reducing dimensionality of data');
[~, ~, v] = svd(stackY);
Ytrans = stackY * v(:, 1:numComponentsEEG);
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

% OLS FIT
betaOls = (X'*X + alphaVal*eye(size(X,2))) \ (X'*Ytrans);

%% FITTING
% disp(['Findig Group-LASSO solution for ' num2str(nLambda) ' regularization (lambda) values']) ;

%Find the optimal regularization parameter for the group-lasso solutions
betaInit = zeros(size(X,2), numComponentsEEG);
betaVal = cell(1, nLambda);
objValues = cell(1, nLambda);
gcvError = zeros(1, nLambda);
df = zeros(1, nLambda);
indexer = cell2mat(arrayfun(@(x) return_index(roiSizes, roiIdx, x), 1:numSubs,'uni',0));
for iLamb = 1:nLambda
    [betaVal{iLamb}, objValues{iLamb}, res] = get_solution_frobenius(X, Ytrans, betaInit, lambdaGrid(iLamb), alphaVal, tol, MAX_ITER, penalties, indices);
    betaInit = betaVal{iLamb};
    betaVal{iLamb} = V * betaVal{iLamb} * v(:,1:numComponentsEEG)'; %transform back to original space (permuted forward matrices)
    rss = norm(stackY-stackedForwards*betaVal{iLamb}(indexer, :), 'fro')^2 / n;
    [gcvError(iLamb), df(iLamb)] = compute_gcv(rss, betaInit, betaOls, grpSizes, n);
end
[~, bestIndex] = min(gcvError);
YhatLASSO = stackedForwards*betaVal{bestIndex}(indexer, :);


%% compute average of metrics
for iSub = 1:numSubs
    range = cell2mat(arrayfun(@(x)  numel(roiIdx{iSub}{x}),1:numROIs,'uni',false));
    range = [0 cumsum(range)];
    temp = betaVal{bestIndex}(return_index(roiSizes, roiIdx, iSub), :);
    roiActivity(:,:,iSub) = cell2mat(arrayfun(@(x) mean(temp(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
end
roiActivityMean(bb,:,:) = mean(roiActivity,3);
roiActivitySum(bb,:,:) = sum(roiActivity,3);
sampleN(bb,:)=pickN;

% %% min-norm
% [betaMinNorm, ~, lambdaMinNorm] = minimum_norm(stackedForwards, Ylo, nLambdaRidge);
% for iSub = 1:numSubs
%     range = cell2mat(arrayfun(@(x)  numel(roiIdx{iSub}{x}),1:numROIs,'uni',false));
%     range = [0 cumsum(range)];
%     tempMinNorm = betaMinNorm(ridgeRange(iSub)+1:ridgeRange(s+1), :);
%     regionActivityMinNorm(:,:,iSub) = cell2mat(arrayfun(@(x) mean(tempMinNorm(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
% end

end % end bootstrap

save('simulOutput2/realDataLasso.mat','sampleN','roiActivityMean','roiActivitySum')

%% Plot ROI activity over time for multiple repeats with CI

% plot different outputs in separate line
count = 1;
figure;set(gcf,'position',[10,10,2400,200])
color = {'r','b'};
% get the 95% CI
ci95 = prctile(roiActivityMean,[2.5 97.5]);
% value for normalising the signal
normVal = max(abs(ci95(:)));
% plot
for iRoi = 1:numROIs
    subplot(1,9,count);hold on
    plot(t, squeeze(mean(roiActivityMean(:,iRoi,:))) / normVal ,color{mod(iRoi,2)+1},'LineWidth',2);
    patch([t fliplr(t)], [squeeze(ci95(1,iRoi,:))' fliplr(squeeze(ci95(2,iRoi,:))')]/ normVal,...
        color{mod(iRoi,2)+1},'FaceAlpha',0.2, 'EdgeColor','none');
    line(t, zeros(size(t)),'Color','k')
    test0 = squeeze(ci95(1,iRoi,:)<=0 & 0<=ci95(2,iRoi,:) ); % test if 0 included in 95CI
    sig0 = find((test0==0)); % get sig indexes
    if mod(iRoi,2)
        scatter(t(sig0), repmat(-0.8,length(sig0),1),'b.')
    else
        scatter(t(sig0), repmat(-0.9,length(sig0),1),'r.')
    end
    tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
    ylim([-1 1]);
    if mod(iRoi,2) == 0
        count=count+1;
    end
end
title(['Lasso' tt(1:end-2)])
legend('left')

saveas(gcf,['figures/realDataLasso'],'png')
saveas(gcf,['figures/realDataLasso' ],'fig')
print(gcf,['figures/realDataLasso' ],'-depsc')

