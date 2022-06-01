clearvars;close all;
% Retrieve sources using real data from Lim 2017

addpath(genpath([pwd filesep 'subfunctions']))
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
% load('/Users/marleneponcet/Documents/Git/templateMinNorm/createTemplate/customEGI128.mat')
% avMap = customTemplates.data;
% listROIs = customTemplates.ROInames;
numROIs = length(listROIs);
nLambdaRidge = 50; 
lowPassNF1 = 1; % filter 1 or not 0
numFq2keep = 5; % nb of harmonics to keep in the signal
totBoot = 0; %500 or 0 for no bootstap

%%%%%%
%% LOAD AND FILTER THE EEG DATA (as specified) 
%%%%%%
sbjList = dir('realData/eegdata/skeri*');
numSubs = length(sbjList);

fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
for iSub=1:numSubs
    clear fwdMatrix roiInfo Axx
    % fwd file
    load(['realData/forwardSolutions/forwardAndRois-' sbjList(iSub).name '.mat'])
    fullFwd{iSub} = fwdMatrix;
    % go through each ROI and save the corresponding fwdMesh values
    % corresponding to the indexes of that ROI
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
        % save the index for each ROI
        idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
    end
    % eegFile
    % PlosOne uses condNmbr = 8;
    Axx = load(['realData/eegdata/' sbjList(iSub).name '/Exp_MATL_HCN_128_Avg/Axx_c008.mat']);
    
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
figure;plot(t,squeeze(mean(Y))','k')
saveas(gcf,['figures/realData' num2str(numFq2keep)],'png')

%%% Re-ref to average
for iSub=1:numSubs
    Y(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end



if totBoot>0
    betaAverage = zeros(totBoot,numROIs,size(Y,3));
    betaAveragePCA = betaAverage;
    retrievePlos = betaAverage;
    retrievePlosMean = betaAverage;
    retrieveWhole = betaAverage;
    retrieveWholeMean = betaAverage;
    
    % bootstrap loop for 95% CI
    
    for bb=1:totBoot
        clear region regionMean regionWhole
        if mod(bb,10)==0;fprintf('boot nb %d',bb);end
        pickN = randi(9,1,9);
        sampleY = Y(pickN,:,:);
        %%%%%% "PCA"
        tempY = permute(sampleY,[2 1 3]);
        stackY = reshape(tempY,[size(Y,2)*numSubs,size(Y,3)]);
        stackY = bsxfun(@minus,stackY, mean(stackY));
        [u1, s1, v1] = svd(stackY);
        numComponents = 5;
        Ylo = u1(:,1:numComponents)*s1(1:numComponents,1:numComponents)*v1(:, 1:numComponents)';
        Y2 = reshape(Ylo,[size(Y,2),numSubs,size(Y,3)]);
        Ypca = permute(Y2,[2 1 3]);
        
     
        %% compute minimum norm
        
        % min_norm on average data: get beta values for each ROI over time
        [betaAverage(bb,:,:), lambda] = minNormFast_lcurve(avMap, squeeze(mean(sampleY,1)));
        [betaAveragePCA(bb,:,:), lambdaPCA] = minNormFast_lcurve(avMap, squeeze(mean(Ypca,1)));
        
        %%%%%%%%%% PlosOne procedure
        stackedFwd=[];stackedFullFwd = [];
        for iSub=pickN
            stackedFwd = blkdiag(stackedFwd, [roiFwd{iSub,:}]);
            stackedFullFwd = blkdiag(stackedFullFwd, [fullFwd{iSub}]);
        end
        stackedFwd = bsxfun(@minus,stackedFwd, mean(stackedFwd));
        stackedFullFwd = bsxfun(@minus,stackedFullFwd, mean(stackedFullFwd));
        [betaPlos,lambdaPlos] = minNormFast(stackedFwd,Ylo,50);
        [betaWhole,lambdaWhole] = minNormFast(stackedFullFwd, Ylo,50);
        
        prevRange = 0; % for counting from prev sbj
        nbS=1; % sbj count
        for iSub=pickN
            % get the number of indexes per ROI for this subj
            rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
            range = [0 cumsum(rangeROI)] + prevRange;
            region(nbS,:,:) = cell2mat(arrayfun(@(x) sum(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
            regionMean(nbS,:,:) = cell2mat(arrayfun(@(x) mean(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
            % whole brain use idxROIfwd & add 20484 for each sbj to account
            % for stacked data
            regionWhole(nbS,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x} + (nbS-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
            regionWholeMean(nbS,:,:) = cell2mat(arrayfun(@(x) mean(betaWhole(idxROIfwd{iSub,x} + (nbS-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
            % increment
            prevRange = range(end);nbS = nbS+1;
        end
        retrievePlos(bb,:,:) = mean(region); % for comparison with sum
        retrievePlosMean(bb,:,:) = mean(regionMean); % = code
        retrieveWhole(bb,:,:) = mean(regionWhole);
        retrieveWholeMean(bb,:,:) = mean(regionWholeMean);
        sampleN(bb,:) = pickN;    
    end
    save('simulOutput2/realData500.mat','sampleN','retrievePlos','retrievePlosMean','betaAverage','betaAveragePCA','retrieveWhole','retrieveWholeMean')
else
    
    %%%%%% "PCA"
    tempY = permute(Y,[2 1 3]);
    stackY = reshape(tempY,[size(Y,2)*numSubs,size(Y,3)]);
    stackY = bsxfun(@minus,stackY, mean(stackY));
    [u1, s1, v1] = svd(stackY);
    numComponents = 5;
    Ylo = u1(:,1:numComponents)*s1(1:numComponents,1:numComponents)*v1(:, 1:numComponents)';
    Y2 = reshape(Ylo,[size(Y,2),numSubs,size(Y,3)]);
    Ypca = permute(Y2,[2 1 3]);
    
    %% compute minimum norm
    
    % min_norm on average data: get beta values for each ROI over time
    [betaAverage, lambda] = minNormFast_lcurve(avMap, squeeze(mean(Y,1)));
    [betaAveragePCA, lambdaPCA] = minNormFast_lcurve(avMap, squeeze(mean(Ypca,1)));

    %%%%%%%%%% PlosOne procedure
    stackedFwd=[]; stackedFullFwd =[];
    for iSub=1:numSubs
        stackedFwd = blkdiag(stackedFwd, [roiFwd{iSub,:}]);
        stackedFullFwd = blkdiag(stackedFullFwd, [fullFwd{iSub}]);
    end
    stackedFwd = bsxfun(@minus,stackedFwd, mean(stackedFwd));
    [betaPlos,lambdaPlos] = minNormFast(stackedFwd,Ylo,50);
    [betaWhole,lambdaWhole] = minNormFast(stackedFullFwd, Ylo,50);
   
    prevRange = 0; % for counting from prev sbj
    for iSub=1:numSubs
        % get the number of indexes per ROI for this subj
        rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
        range = [0 cumsum(rangeROI)] + prevRange;
        region(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
        regionMean(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
        prevRange = range(end);
        % whole brain use idxROIfwd & add 20484 for each sbj to account
        % for stacked data
        regionWhole(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaWhole(idxROIfwd{iSub,x} + (iSub-1)*(length(betaWhole)/numSubs),:)),1:numROIs,'uni',false)');
    end
    retrievePlos = mean(region); % for comparison with sum
    retrievePlosMean = mean(regionMean); % = code
    retrieveWhole = mean(regionWhole);
end

% save('simulOutput/realDataCustom.mat','retrievePlos','retrievePlosMean','betaAverage','betaAveragePCA','betaMinNorm')


%%%%% figures ROIs
if totBoot>0 
    % plot different outputs in separate line
    count = 1; 
    plotMod = 1; figure;set(gcf,'position',[10,10,2400,600])
%     plotMod = 4; figure;set(gcf,'position',[10,10,2400,1000])
    color = {'r','b'};
    for oo=1:plotMod
        if oo==1
            data = betaAveragePCA;tname = 'Template';
        elseif oo==2
            data = retrieveWholeMean;tname = 'Individual';
        elseif oo==3
            data = retrievePlosMean;tname = 'Individual-ROI';
%             data = betaAverage1;tname = 'template';
        elseif oo==4
            data = retrievePlos;tname = 'constrainedSum';
        end
        % get the 95% CI
        ci95 = prctile(data,[2.5 97.5]);
        % value for normalising the signal
        normVal = max(abs(ci95(:)));
        % plot
        for iRoi = 1:numROIs
            subplot(plotMod,9,count);hold on
            plot(t, squeeze(mean(data(:,iRoi,:))) / normVal ,color{mod(iRoi,2)+1},'LineWidth',2);
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
        title([tname tt(1:end-2)])
    end
    legend('left')
    saveas(gcf,['figures/realDataCustom'],'png')
    saveas(gcf,['figures/realDataCustom' ],'fig')
    print(gcf,['figures/realDataCustom' ],'-depsc')   
    
else
    %% Plot without bootstrap
    count = 1;
    figure;set(gcf,'position',[100,100,2200,700])
    % template
    for iRoi = 1:2:numROIs
        % need to normalise the signal
        subplot(2,9,count);hold on
        plot(t, betaAverage(iRoi,:) / max(max(abs(betaAverage))) ,'LineWidth',2);
        plot(t, betaAverage(iRoi+1,:) / max(abs(betaAverage(:))) ,'LineWidth',2);
        line(t, zeros(size(t)),'Color','k')
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    title(['template' tt(1:end-2)])
    for iRoi = 1:2:numROIs
        % need to normalise the signal
        subplot(2,9,count);hold on
        plot(t, retrievePlosMean(iRoi,:) / max(max(abs(retrievePlosMean))) ,'LineWidth',2);
        plot(t, retrievePlosMean(iRoi+1,:) / max(max(abs(retrievePlosMean))) ,'LineWidth',2);
        line(t, zeros(size(t)),'Color','k')
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
    title(['constrained' tt(1:end-2)])
    legend('left','right')
    saveas(gcf,'figures/realData','png')
    saveas(gcf,'figures/realData','fig')
    print(gcf,'figures/realData','-depsc')
    
    % template PCA
    for iRoi = 1:2:numROIs
        % need to normalise the signal
        subplot(3,9,count);hold on
        plot(t, betaAveragePCA(iRoi,:) / max(max(abs(betaAveragePCA))) ,'LineWidth',2);
        plot(t, betaAveragePCA(iRoi+1,:) / max(abs(betaAveragePCA(:))) ,'LineWidth',2);
        line(t, zeros(size(t)),'Color','k')
        tt = cell2mat(listROIs(iRoi));title(tt(1:end-2))
        ylim([-1 1]);count=count+1;
    end
end


if totBoot == 0
addpath('/Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated/')
addpath(genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork'));

load('simulOutput2/dataLasso.mat')
reLasso = squeeze(mean(unstackedYhatLASSO,2));
%%% topography 250 ms (as in Plos One)
reBeta = avMap * betaAveragePCA; % reconstruct elec amplitude * time
% for individual source loc, need to reconstruct per sbj
YhatPlos=stackedFwd*betaPlos;
unstackedYhatPlos = reshape(YhatPlos,128,numSubs,780);
rePlos = squeeze(mean(unstackedYhatPlos,2));
YhatWhole=stackedFullFwd*betaWhole;
unstackedYhatWhole = reshape(YhatWhole,128,numSubs,780);
reWhole = squeeze(mean(unstackedYhatWhole,2));

timeToPlot = 250;
layoutData = '/Users/marleneponcet/Documents/Git/templateMinNorm/createTemplate/layout/GSN-HydroCel-128.sfp';
mm = round(max([abs(mean(Ypca(:,:,timeToPlot),1))'; abs(reBeta(:,timeToPlot)); ...
    abs(rePlos(:,timeToPlot)); abs(reLasso(:,timeToPlot))]),1);

% 2D fieldtrip
figure;
subplot(2,3,1);plotTopo(mean(Ypca(:,:,timeToPlot),1),layoutData);caxis([-mm mm]);title('data')
subplot(2,3,2);plotTopo(reBeta(:,timeToPlot),layoutData);caxis([-mm mm]);title('Template')
subplot(2,3,3);plotTopo(rePlos(:,timeToPlot),layoutData);caxis([-mm mm]);title('Individual')
subplot(2,3,4);plotTopo(reWhole(:,timeToPlot),layoutData);caxis([-mm mm]);title('Individual Subset')
subplot(2,3,5);plotTopo(reWhole(:,timeToPlot),layoutData);caxis([-mm mm]);title('Lasso')
colorcet('D1')
set(gcf, 'Position', [100 500 1000 600]);
saveas(gcf,['figures/realDataTopo2D'],'png')
saveas(gcf,['figures/realDataTopo2D'],'fig')
print(gcf,['figures/realDataTopo2D'],'-depsc')

% % 2D Justin
% figure;
% subplot(2,2,1);plotOnEgi(mean(Ypca(:,:,timeToPlot),1));caxis([-mm mm]);title('data')
% subplot(2,2,2);plotOnEgi(reBeta(:,timeToPlot));caxis([-mm mm]);title('template')
% subplot(2,2,3);plotOnEgi(rePlos(:,timeToPlot));caxis([-mm mm]);title('constrained')
% subplot(2,2,4);plotOnEgi(reWhole(:,timeToPlot));caxis([-mm mm]);title('tailored')
% colorcet('D1')
% set(gcf, 'Position', [100 500 1000 600]);
% saveas(gcf,['figures/realDataTopo2D_justin'],'png')
% saveas(gcf,['figures/realDataTopo2D_justin'],'fig')
% print(gcf,['figures/realDataTopo2D_justin'],'-depsc')

% 3D
figure;
subplot(2,3,1);
plotContourOnScalp(mean(Ypca(:,:,timeToPlot),1),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('data')
subplot(2,3,2);
plotContourOnScalp(reBeta(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('template')
subplot(2,3,3);
plotContourOnScalp(rePlos(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('constrained')
subplot(2,3,4);
plotContourOnScalp(reWhole(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('tailored')
subplot(2,3,5);
plotContourOnScalp(reWhole(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
set(gcf, 'Position', [100 500 1000 600]);
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('lasso')
saveas(gcf,['figures/realDataTopo3D'],'png')
saveas(gcf,['figures/realDataTopo3D'],'fig')
print(gcf,['figures/realDataTopo3D'],'-depsc')
end

%%% RMSE
% diff3 = squeeze(mean(Ypca(:,1,:),1)) - test(1,:)';
% elecRMSE = rms(squeeze(mean(Ypca(:,1,:),1)) - test(1,:)');
% diff = squeeze(mean(Ypca,1)) - test;
betaRMSEpca = rms((squeeze(mean(Ypca,1)) - reBeta));
plosRMSEpca = rms((squeeze(mean(Ypca,1)) - rePlos));
betaRMSE = rms((squeeze(mean(Y,1)) - reBeta));
plosRMSE = rms((squeeze(mean(Y,1)) - rePlos));
figure;hold on
plot(t,betaRMSEpca,'LineWidth',2); plot(t,plosRMSEpca,'LineWidth',2); 
plot(t,betaRMSE,'LineWidth',2); plot(t,plosRMSE,'LineWidth',2);
legend('tempYpca','plosYpca','tempY','plosY')
xlabel('time');ylabel('RMSE');
saveas(gcf,'figures/realDataRMSE','png')
