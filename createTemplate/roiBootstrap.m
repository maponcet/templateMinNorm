%%%%% ROI template bootstrap with replacement for increasing nb of max sbj
clearvars;
listFiles = dir('/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/forward*');
load('averageMap50.mat')
nbBoot = 1000;
rmseROI = zeros(nbBoot,length(listFiles),size(avMap.activity,2));

for maxSbj = 1:length(listFiles)
    for bb=1:nbBoot
        if mod(bb,100)==0
            fprintf('maxSbj%d boostrap%d \n' ,maxSbj,bb)
        end
        
        % pick randomly maxSbj sbj with replacement
        pickSS = randi(length(listFiles),1,maxSbj);
        roiMap = zeros(size(avMap.activity,1),size(avMap.activity,2),length(pickSS));
        
        for iSubj = 1:length(pickSS)
            clear fwdMatrix roiInfo
            load([listFiles(pickSS(iSubj)).folder filesep listFiles(pickSS(iSubj)).name])
            % go through each ROI and average all the mesh indexes corresponding to that ROI
            for rr=1:length(avMap.roiNames)
                clear indexROI
                indexROI = find(strcmp(avMap.roiNames(rr),{roiInfo.name}));
                roiMap(:,rr,iSubj) = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
            end
            
        end
        rmseROI(bb,maxSbj,:) =  rms(mean(roiMap,3) - avMap.activity);
    end
end
listROIs = avMap.roiNames;
save roiBootstrapFromfwd.mat rmseROI listROIs

%%%%%%%%%%%%%%%%%%%%%%%
% rmse: boot, sbj, roi
% compute %error: rmse/rms(avMap.activity)
load('roiBootstrapFromfwd.mat')
load('averageMap50.mat')
prctError=zeros(size(rmseROI));
for bb=1:size(rmseROI,1)
    prctError(bb,:,:) = squeeze(rmseROI(bb,:,:))./rms(avMap.activity)*100;
end
% plot rmse/rms

figure;hold on;
set(gcf,'position',[100,100,1500,700])
for src=1:size(prctError,3)
    bb = prctile(prctError(:,:,src),97.5);
    subplot(3,6,src)
    plot(mean(prctError(:,:,src),1),'LineWidth',2)
    patch([1:size(prctError,2) size(prctError,2):-1:1],[ prctile(prctError(:,:,src),2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
    xlabel('nb sbj included')
    ylabel('% error')
    set(gca, 'YScale', 'log','XScale', 'log')
    xlim([0 50]);ylim([0 250])
    [aa, bb] = polyfit(log10(1:50),log10(mean(prctError(:,:,src),1)),1);
    title([listROIs{src} ' S=' num2str(aa(1),'%.2f')])
end
saveas(gcf,['figures' filesep 'prctErrorLog'],'png')
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

