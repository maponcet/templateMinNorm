%%%%% Percent error in ROI compared to template for at each N for increasing N

clearvars;
listFiles = dir('/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/forward*');
load('templates/template_EGI128.mat')
avMap = templates.weights;
listROIs = templates.listROIs;
nbBoot = 800;
percErr = zeros(nbBoot,length(listFiles),size(avMap,2));

for maxSbj = 1:length(listFiles)
    for bb=1:nbBoot
        if mod(bb,10)==0
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
        percErr(bb,maxSbj,:) =  arrayfun(@(roi) (mean(roiMap(:,roi,:),3) - avMap(:,roi))'.^2 / avMap(:,roi)'.^2, 1:length(listROIs));        
    end
end
save roiPrcErr.mat percErr listROIs

%%%%%%%%%%%%%%%%%%%%%%%
% percErr: boot, sbj, electrodes, roi
load('roiPrcErr.mat')

gg=1;
figure;hold on;
set(gcf,'position',[100,100,1800,500])
for src=[1:2:18 2:2:18]
    bb = prctile(percErr(:,:,src),97.5);
    subplot(2,9,gg); hold on;
%     plot(1./(1:size(percErr,2)),'LineWidth',2)
%     plot(mean(sqrt(percErr(:,:,src)),1),'LineWidth',2)
%     plot(1./sqrt(1:size(percErr,2)),'LineWidth',2)
    plot(mean(percErr(:,:,src),1),'LineWidth',2)
    patch([1:size(percErr,2) size(percErr,2):-1:1],[ prctile(percErr(:,:,src),2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
    xlabel('log (N)')
    ylabel('log (% error)')
    set(gca, 'YScale', 'log','XScale', 'log')
    xlim([1 50]); ylim([10^-3 10])
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


