% create average activity map from the forward .mat files from the 25 sbj
% in the PlosOne paper + the 50 total sbj (average ROI locations)
addpath('/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/requiredFunctions/');
addpath(genpath('/Volumes/Amrutam/Marlene/Git/svndl_code/'));
setpref('mrLASSO','scalpFileDir','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/anatomy');



%%% for the sbj in PlosOne paper
listFiles = dir('/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/forwardSolutions/forward*');

% plosOne: we defined the detailed 3D shape of each of 18 visual ROIs in 25 participants 
% (V1-L, V1-R, V2v-L, V2v-R, V2d-L, V2d-R, V3v-L, V3v-R, V3d-L, V3d-R, ...
% V4-L, V4-R, V3A-L, V3A-R, LOC-L, LOC-R, MT-L, MT-R)
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};

% LOC missing for 60, 72 (after removing them = 25 participants which is
% the same as in the paper)
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

% average the roiMap accross sbj
avMap.activity = mean(roiMap,3);
avMap.roiNames = listROIs;

figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotContourOnScalp(avMap.activity(:,roi),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
end
saveas(gcf,'averageMap25Plos','png')
save('averageMap25Plos.mat','avMap')





%%%% For the 50 sbj
clearvars
listFiles = dir('/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/forward*');
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

% average the roiMap accross sbj
% avMap.activity = mean(roiMap,3);
avMap = mean(roiMapSum,3);
figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotContourOnScalp(avMap(:,roi),'skeri0044','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
end
saveas(gcf,'averageMap50','png')
save('averageMap50Sum.mat','avMap','listROIs')
avMap= mean(roiMapMean,3);
save('averageMap50Mean.mat','avMap','listROIs')

figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMap.sum(:,roi))
    caxis([-8900 8900])
%     colorbar
end
saveas(gcf,'figures/averageMap50sameRangeSUM','png')

figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMap.mean(:,roi))
    caxis([-150 150])
end
saveas(gcf,'figures/averageMap50sameRangeMEAN','png')





%Scale cdata -1 -> 1
cdata =  (cdata/max(abs(cdata)));
%Push -1:1 range to 0:1
cdata = (cdata+1)/2;

figure('position', [200, 0, 1500, 800])
aa=1;
for roi=1:2:18
    subplot(3,3,aa)
    title(listROIs(roi))
    plotOnEgi(avMap.activity(:,roi) ./ max(abs(avMap.activity(:,1))))
    colorbar;
    aa=aa+1;
end 
saveas(gcf,'averageMap50scaleV1','png')
figure('position', [200, 0, 1500, 800])
aa=1;
for roi=1:2:18
    subplot(3,3,aa)
    title(listROIs(roi))
    plotOnEgi(avMap.activity(:,roi) ./ max(abs(avMap.activity(:,1))))
    caxis([-2 2]);colorbar;
    aa=aa+1;
end
saveas(gcf,'averageMap50scaleV1sameRange','png')

% % NONSENSE!
% figure;plotOnEgi(avMap.activity(:,1))
% cmap = jmaColors('arizona');
% caxis([-50 150])
% % 32 colours from 0 to 150
% % remove 20 colours in the first 32 cmap
% cmap2 = cmap([1:3:32 32:end],:);
% colormap(cmap2);
% colorbar
%   
% figure('position', [200, 0, 1500, 800])
% for roi=1:18
%     subplot(3,6,roi)
%     title(listROIs(roi))
%     plotOnEgi(avMap.activity(:,roi))
%     caxis([-50 150])
%     colormap(cmap2);
% end
