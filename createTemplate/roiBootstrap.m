%%%%% ROI template bootstrap with replacement for increasing nb of max sbj

clearvars;
listFiles = dir('/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/forward*');
load('averageMap50Sum.mat')
nbBoot = 20;
percErr = zeros(nbBoot,length(listFiles),size(avMap,1),size(avMap,2));

for maxSbj = 1:length(listFiles)
    for bb=1:nbBoot
        if mod(bb,100)==0
            fprintf('maxSbj%d boostrap%d \n' ,maxSbj,bb)
        end
        
        % pick randomly sbj with replacement
        pickSS = randi(length(listFiles),1,maxSbj);
        roiMap = zeros(size(avMap,1),size(avMap,2),length(pickSS));
        
        for iSubj = 1:length(pickSS)
            clear fwdMatrix roiInfo
            load([listFiles(pickSS(iSubj)).folder filesep listFiles(pickSS(iSubj)).name])
            % go through each ROI and average all the mesh indexes corresponding to that ROI
            for rr=1:length(listROIs)
                clear indexROI
                indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
                roiMap(:,rr,iSubj) = sum(fwdMatrix(:,roiInfo(indexROI).meshIndices),2); 
            end
            
        end
        percErr(bb,maxSbj,:,:) =  (mean(roiMap,3) - avMap).^2 ./ avMap.^2 ;
    end
end
save roiBootstrapFromfwd.mat percErr listROIs

%%%%%%%%%%%%%%%%%%%%%%%
% percErr: boot, sbj, electrodes, roi
load('roiBootstrapFromfwd.mat')
prctError = rmseROI;
% % average across electrodes
% prctError = squeeze(mean(percErr,3));

gg=1;
figure;hold on;
set(gcf,'position',[100,100,1800,500])
for src=[1:2:18 2:2:18]
    bb = prctile(prctError(:,:,src),97.5);
    subplot(2,9,gg); hold on;
    plot(mean(prctError(:,:,src),1),'LineWidth',2)
%     plot(sqrt(mean(prctError(:,:,src),1)),'LineWidth',2)
    patch([1:size(prctError,2) size(prctError,2):-1:1],[ prctile(prctError(:,:,src),2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
    xlabel('log (N)')
    ylabel('log (% error)')
    set(gca, 'YScale', 'log','XScale', 'log')
    xlim([1 50]); ylim([1 60])
%     [aa, bb] = polyfit(log10(1:50),log10(mean(prctError(:,:,src),1)),1);
%     title([listROIs{src} ' S=' num2str(aa(1),'%.2f')])
    title(listROIs{src})
    gg=gg+1;
end
saveas(gcf,['figures' filesep 'prctErrorLog'],'png')
saveas(gcf,['figures' filesep 'prctErrorLog'],'fig')
print(gcf,['figures' filesep 'prctErrorLog'],'-depsc')
% slope = -0.5 so for 1 Yunit needs 2Xunit, if log10 then change 10^1 error needs
% 10^2 sbj; for decrease of half error (log2^1) needs 2^2 = quadruple sbj
% is it just because of the way we do the bootstrap???????


%%%% plot std
figure;hold on;
set(gcf,'position',[100,100,1500,700])
for src=1:size(prctError,3)
    subplot(3,6,src)
    plot(std(prctError(:,:,src)))
    xlabel('nb sbj included')
    ylabel('std')
    xlim([0 50]);ylim([0 50])
%     set(gca, 'YScale', 'log','XScale', 'log')
    [aa, bb] = polyfit(log(1:50),log(std(prctError(:,:,src))),1);
    title([listROIs{src} ' S=' num2str(aa(1),'%.2f')])
end
saveas(gcf,['figures' filesep 'prctErrorStdLog'],'png')

%%%% plot change
figure;hold on;
set(gcf,'position',[100,100,1500,900])
for src=1:size(prctError,3)
    change = prctError(:,1:49,src) - prctError(:,2:50,src);
    bb = prctile(change,97.5);
    subplot(3,6,src)
    plot(mean(change,1),'LineWidth',2)
    title(listROIs{src})
    xlabel('nb sbj included')
    ylabel('change in prctError')
%     line([0 50],[0.001 0.001],'Color','r')
%     ylim([0 2])
    patch([1:size(change,2) size(change,2):-1:1],[ prctile(change,2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
%     ylim([0 0.25])
    xlim([0 50])
end
saveas(gcf,['figures' filesep 'prctErrorChange'],'png')

%%%% plot variability
figure;hold on;
set(gcf,'position',[100,100,1500,700])
for src=1:size(prctError,3)
    bb = prctile(prctError(:,:,src),97.5) - prctile(prctError(:,:,src),2.5);
    subplot(3,6,src)
    plot(mean(prctError(:,:,src),1),'LineWidth',2)
    title(listROIs{src})
    xlabel('nb sbj included')
    ylabel('97.5-2.5 prctError prctile')
    ylim([0 1])
    xlim([0 50])
end
saveas(gcf,['figures' filesep 'prctErrorVar'],'png')

%%% previous plots
%%%%%%%%%%%%%%%%%%%%%%%
%% Plot average-based ROIs RMSE (+variability)
figure;hold on;
set(gcf,'position',[100,100,1500,700])
for src=1:size(rmseROI,3)
    bb = prctile(rmseROI(:,:,src),97.5);
    subplot(3,6,src)
    plot(mean(rmseROI(:,:,src),1),'LineWidth',2)
    patch([1:size(rmseROI,2) size(rmseROI,2):-1:1],[ prctile(rmseROI(:,:,src),2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
    title(listROIs{src})
    xlabel('nb sbj included')
    ylabel('rmse')
%     ylim([0 0.3])
    xlim([0 50])
end
saveas(gcf,['figures' filesep 'rmseAvRoiv2'],'png')

%%%% plot change
figure;hold on;
set(gcf,'position',[100,100,1500,900])
for src=1:size(rmseROI,3)
    change = rmseROI(:,1:49,src) - rmseROI(:,2:50,src);
    bb = prctile(change,97.5);
    subplot(3,6,src)
    plot(mean(change,1),'LineWidth',2)
    title(listROIs{src})
    xlabel('nb sbj included')
    ylabel('change in RMSE')
%     line([0 50],[0.001 0.001],'Color','r')
%     ylim([0 0.045])
    patch([1:size(change,2) size(change,2):-1:1],[ prctile(change,2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
%     ylim([0 0.25])
    xlim([0 50])
end
saveas(gcf,['figures' filesep 'rmseChange'],'png')

%%% plot variability 
figure;hold on;
set(gcf,'position',[100,100,1500,700])
for src=1:size(rmseROI,3)
    bb = prctile(rmseROI(:,:,src),97.5) - prctile(rmseROI(:,:,src),2.5);
    subplot(3,6,src)
    plot(mean(rmseROI(:,:,src),1),'LineWidth',2)
    title(listROIs{src})
    xlabel('nb sbj included')
    ylabel('97.5-2.5 rmse prctile')
%     ylim([0 0.15])
    xlim([0 50])
end
saveas(gcf,['figures' filesep 'rmseVar'],'png')

