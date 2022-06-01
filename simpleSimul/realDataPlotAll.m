clearvars

load('simulOutput2/realData100.mat')
betaAverageAll = betaAverage;
retrievePlosMeanAll = retrievePlosMean;
retrieveWholeMeanAll = retrieveWholeMean;
retrievePlosAll = retrievePlos;
retrieveWholeAll = retrieveWhole;

load('simulOutput2/realData200.mat')
betaAverageAll = [betaAverageAll; betaAveragePCA];
retrievePlosMeanAll = [retrievePlosMeanAll; retrievePlosMean];
retrieveWholeMeanAll = [retrieveWholeMeanAll; retrieveWholeMean];
retrievePlosAll = [retrievePlosAll; retrievePlos];
retrieveWholeAll = [retrieveWholeAll; retrieveWhole];

load('simulOutput2/realData300.mat')
betaAverageAll = [betaAverageAll; betaAveragePCA];
retrievePlosMeanAll = [retrievePlosMeanAll; retrievePlosMean];
retrieveWholeMeanAll = [retrieveWholeMeanAll; retrieveWholeMean];
retrievePlosAll = [retrievePlosAll; retrievePlos];
retrieveWholeAll = [retrieveWholeAll; retrieveWhole];

load('simulOutput2/realData400.mat')
betaAverageAll = [betaAverageAll; betaAveragePCA];
retrievePlosMeanAll = [retrievePlosMeanAll; retrievePlosMean];
retrieveWholeMeanAll = [retrieveWholeMeanAll; retrieveWholeMean];
retrievePlosAll = [retrievePlosAll; retrievePlos];
retrieveWholeAll = [retrieveWholeAll; retrieveWhole];

load('simulOutput2/realData500.mat')
betaAverageAll = [betaAverageAll; betaAveragePCA];
retrievePlosMeanAll = [retrievePlosMeanAll; retrievePlosMean];
retrieveWholeMeanAll = [retrieveWholeMeanAll; retrieveWholeMean];
retrievePlosAll = [retrievePlosAll; retrievePlos];
retrieveWholeAll = [retrieveWholeAll; retrieveWhole];

load('simulOutput2/realDataLasso.mat')
% plot different outputs in separate line
load('simulOutput2/timeline.mat')
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
count = 1;totBoot = 500;plotMod=4;
figure;set(gcf,'position',[10,10,2400,1000])
color = {'r','b'};
for oo=1:plotMod
    if oo==1
        data = betaAverageAll;tname = 'Template';
    elseif oo==2
        data = retrieveWholeMeanAll;tname = 'Individual';
    elseif oo==3
        data = retrievePlosMeanAll;tname = 'Individual-Subset';
    elseif oo==4
        data = roiActivityMean;tname = 'Lasso';
    elseif oo==5
        data = retrieveWholeAll;tname = 'IndividualSum';
    elseif oo==6
        data = retrievePlosAll;tname = 'IndividualSum-Subset';
    end
    % get the 95% CI
    ci95 = prctile(data,[2.5 97.5]);
    % value for normalising the signal
    normVal = max(abs(ci95(:)));
    % plot
    for iRoi = 1:length(listROIs)
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
saveas(gcf,['figures/realDataAll'],'png')
saveas(gcf,['figures/realDataAll' ],'fig')
print(gcf,['figures/realDataAll' ],'-depsc')