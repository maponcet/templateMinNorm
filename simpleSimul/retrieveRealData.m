clearvars;close all;
% Retrieve sources using real data from Lim 2017

addpath(genpath([pwd filesep 'subfunctions']))
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
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
    
    % bootstrap loop for 95% CI
    
    for bb=1:totBoot
        clear region regionMean
        
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
        stackedForwards=[];
        for iSub=pickN
            stackedForwards = blkdiag(stackedForwards, [roiFwd{iSub,:}]);
        end
        stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));
        [betaPlos,lambdaPlos] = minNormFast(stackedForwards,Ylo,50);
        
        prevRange = 0; % for counting from prev sbj
        nbS=1; % sbj count
        for iSub=pickN
            % get the number of indexes per ROI for this subj
            rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
            range = [0 cumsum(rangeROI)] + prevRange;
            region(nbS,:,:) = cell2mat(arrayfun(@(x) sum(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
            regionMean(nbS,:,:) = cell2mat(arrayfun(@(x) mean(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
            prevRange = range(end);nbS = nbS+1;
        end
        retrievePlos(bb,:,:) = mean(region); % for comparison with sum
        retrievePlosMean(bb,:,:) = mean(regionMean); % = code
        
    end
    
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
    stackedForwards=[];
    for iSub=numSubs
        stackedForwards = blkdiag(stackedForwards, [roiFwd{iSub,:}]);
    end
    stackedForwards = bsxfun(@minus,stackedForwards, mean(stackedForwards));
    [betaPlos,lambdaPlos] = minNormFast(stackedForwards,Ylo,50);
    
    prevRange = 0; % for counting from prev sbj
    for iSub=numSubs
        % get the number of indexes per ROI for this subj
        rangeROI = cell2mat(arrayfun(@(x)  numel(idxROIfwd{iSub,x}),1:numROIs,'uni',false));
        range = [0 cumsum(rangeROI)] + prevRange;
        region(iSub,:,:) = cell2mat(arrayfun(@(x) sum(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
        regionMean(iSub,:,:) = cell2mat(arrayfun(@(x) mean(betaPlos(range(x)+1:range(x+1), :)),1:numROIs,'uni',false)');
        prevRange = range(end);
    end
    retrievePlos = mean(region); % for comparison with sum
    retrievePlosMean = mean(regionMean); % = code
end

save('simulOutput/realData.mat','retrievePlos','retrievePlosMean','betaAverage','betaAveragePCA')

reBeta = avMap * squeeze(mean(betaAveragePCA)); % reconstruct elec amplitude * time
rePlos = avMap * squeeze(mean(retrievePlosMean)); % reconstruct elec amplitude * time
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

%%%%% figures
if totBoot>0 
    % plot different outputs in separate line
    count = 1; 
    plotMod = 2; figure;set(gcf,'position',[10,10,2400,600])
%     plotMod = 4; figure;set(gcf,'position',[10,10,2400,1000])
    color = {'r','b'};
    for oo=1:plotMod
        if oo==1
            data = betaAveragePCA;tname = 'templatePCA';
        elseif oo==2
            data = retrievePlosMean;tname = 'constrainedMean';
        elseif oo==3
            data = betaAverage;tname = 'template';
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
            test0 = zeros(length(data),1);
            for tt=1:length(data)
                test0(tt) = sum(data(:,iRoi,tt)>0) / totBoot;
                if test0(tt) > 0.5
                    test0(tt) = 1-test0(tt);
                end
            end
            sig0 = find((test0<0.025));
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
    saveas(gcf,['figures/realDataCI' num2str(plotMod)],'png')
    saveas(gcf,['figures/realDataCI' num2str(plotMod)],'fig')
    print(gcf,['figures/realDataCI' num2str(plotMod)],'-depsc')   
    
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


