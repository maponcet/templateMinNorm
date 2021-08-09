

%%% ATTENTION: should plot topo according to the system layout!!!
clearvars
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork


listFiles = dir('/Users/marleneponcet/Documents/data/skeriDATA/forwardBiosemi128/*.mat');
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
aa=0;
roiMapMean=[];roiMapSum=[];
for ff=1:length(listFiles)
    clear fwdMatrix roiInfo
    load([listFiles(ff).folder filesep listFiles(ff).name])
    % check that all rois are there
    if sum(ismember(listROIs, {roiInfo.name})) ~= length(listROIs)
        fprintf('nb of ROIs mismatch for %s \n',listFiles(ff).name)
    else
        aa = aa+1;
        % go through each ROI and average all the mesh indexes corresponding to that ROI
        for rr=1:length(listROIs)
            clear indexROI
            indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            roiMapMean(:,rr,aa) = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
            roiMapSum(:,rr,aa) = sum(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
        end
    end
end

avMap = mean(roiMapSum,3);
max(max(abs(avMap)))
mm = 9750;
figure('position', [200, 0, 1500, 800])
colormap jet
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotTopo(avMap(:,roi),'biosemi128.lay')
%     caxis([-mm mm])
%     colorbar
end
saveas(gcf,'figures/averageMapBiosemi128','png')
save('averageMapBiosemi128.mat','avMap','listROIs')




listFiles = dir('/Users/marleneponcet/Documents/data/skeriDATA/forwardBiosemi64/*.mat');
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
roiMap=[];aa=0;
for ff=1:length(listFiles)
    clear fwdMatrix roiInfo
    load([listFiles(ff).folder filesep listFiles(ff).name])
    % check that all rois are there
    if sum(ismember(listROIs, {roiInfo.name})) ~= length(listROIs)
        fprintf('nb of ROIs mismatch for %s \n',listFiles(ff).name)
    else
        aa = aa+1;
        % go through each ROI and average all the mesh indexes corresponding to that ROI
        for rr=1:length(listROIs)
            clear indexROI
            indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            roiMap(:,rr,aa) = sum(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
        end
    end
end

avMap = mean(roiMap,3);
max(max(abs(avMap)))
mm=9100;
figure('position', [200, 0, 1500, 800])
colormap jet
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotTopo(avMap(:,roi),'biosemi64.lay')
    caxis([-mm mm])
%     colorbar
end
saveas(gcf,'figures/averageMapBiosemi64','png')
save('averageMapBiosemi64.mat','avMap','listROIs')





listFiles = dir('/Users/marleneponcet/Documents/data/skeriDATA/test/*.mat');
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
roiMap=[];aa=0;
roiMapMean=[];roiMapSum=[];
for ff=1:length(listFiles)
    clear fwdMatrix roiInfo
    load([listFiles(ff).folder filesep listFiles(ff).name])
    % check that all rois are there
    if sum(ismember(listROIs, {roiInfo.name})) ~= length(listROIs)
        fprintf('nb of ROIs mismatch for %s \n',listFiles(ff).name)
    else
        aa = aa+1;
        % go through each ROI and average all the mesh indexes corresponding to that ROI
        for rr=1:length(listROIs)
            clear indexROI
            indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            roiMapMean(:,rr,aa) = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
            roiMapSum(:,rr,aa) = sum(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
        end
    end
end
figure;plotTopo(roiMap(1,1),'biosemi64.lay')
