% create average activity map from the forward .mat files from the 25 sbj
% in the PlosOne paper + the 50 total sbj (average ROI locations)
addpath('/Volumes/Amrutam/Marlene/JUSTIN/avROImap/PlosOne/github-archive/requiredFunctions/');
setpref('mrLASSO','scalpFileDir','/Volumes/Amrutam/Marlene/JUSTIN/avROImap/PlosOne/github-archive/datafiles/anatomy');



%%% for the sbj in PlosOne paper
listFiles = dir('/Volumes/Amrutam/Marlene/JUSTIN/avROImap/PlosOne/github-archive/datafiles/forwardSolutions/forward*');

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
            roiMap(:,rr,aa) = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
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
    plotContourOnScalp(avMap.activity(:,roi),'skeri0044','datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
end
saveas(gcf,['figures/' 'averageMap25'],'png')
save('averageMap25.mat','avMap')





%%%% For the 50 sbj
clearvars
listFiles = dir('/Volumes/Amrutam/Marlene/JUSTIN/avROImap/PlosOne/github-archive/forwardAllEGI/forward*');
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
            roiMap(:,rr,aa) = mean(fwdMatrix(:,roiInfo(indexROI).meshIndices),2);
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
    plotContourOnScalp(avMap.activity(:,roi),'skeri0044','datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
end
saveas(gcf,['figures/' 'averageMap50'],'png')
save('averageMap50.mat','avMap')




