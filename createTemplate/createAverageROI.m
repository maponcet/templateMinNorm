% create average activity map from the forward .mat files from the 25 sbj
% in the PlosOne paper + the 50 total sbj (average ROI locations)
% Plot individual ROIs

% addpath('/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/requiredFunctions/');
addpath(genpath('/Users/marleneponcet/Documents/Git/svndl_code/'));
addpath(genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork'));
% setpref('mrLASSO','scalpFileDir','/Volumes/Amrutam/Marlene/JUSTIN/PlosOne/github-archive/datafiles/anatomy');


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
% listFiles = dir('/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/forward*');
listFiles = dir('/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/*.mat');
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
figure; pict = 1;
for roi=2:2:18
    subplot(1,9,pict)
    title(listROIs(roi))
    plotContourOnScalp(avMap(:,roi),'skeri0044','/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
    pict = pict+1;
end
set(gcf,'position', [200, 0, 2000, 300])
saveas(gcf,'averageMap50_3Dturn','png')

save('averageMap50Sum.mat','avMap','listROIs')
avMap= mean(roiMapMean,3);
save('averageMap50Mean.mat','avMap','listROIs')

figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMap(:,roi))
    caxis([-8900 8900])
%     colorbar
end
saveas(gcf,'figures/averageMap50sameRangeSUM','png')

figure('position', [200, 0, 1500, 800])
for roi=1:18
    subplot(3,6,roi)
    title(listROIs(roi))
    plotOnEgi(avMap(:,roi))
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



% plot each ind V1
pickROI = [1 2 17 18];
for rr=1:length(pickROI)
    minScale = min(min(roiMapSum(:,pickROI(rr),:)));
    maxScale = max(max(roiMapSum(:,pickROI(rr),:)));
    scale = max([abs(minScale) abs(maxScale)]);
    figure('position', [200, 100, 1800, 1200])
    for ff=1:25
        subplot(5,5,ff);
        plotOnEgi(roiMapSum(:,pickROI(rr),ff))
        caxis([-scale scale]);
    end
    saveas(gcf,['indROI/scale' listROIs{pickROI(rr)} '_S1-25'],'png')
    figure('position', [200, 100, 1800, 1200])
    for ff=26:length(listFiles)
        subplot(5,5,ff-25);
        plotOnEgi(roiMapSum(:,pickROI(rr),ff))
        caxis([-scale scale]);
    end
    saveas(gcf,['indROI/scale' listROIs{pickROI(rr)} '_S26-50'],'png')
    figure('position', [200, 100, 1800, 1200])
    for ff=1:25
        subplot(5,5,ff);
        plotOnEgi(roiMapSum(:,pickROI(rr),ff))
        colorbar;
    end
    saveas(gcf,['indROI/unscaled' listROIs{pickROI(rr)} '_S1-25'],'png')
    figure('position', [200, 100, 1800, 1200])
    for ff=26:length(listFiles)
        subplot(5,5,ff-25);
        plotOnEgi(roiMapSum(:,pickROI(rr),ff))
        colorbar;
    end
    saveas(gcf,['indROI/unscale' listROIs{pickROI(rr)} '_S26-50'],'png')
end

% plot different representative response
pickROI = 2; pickInd = [2 30 19];
minScale = min(min(roiMapSum(:,pickROI,pickInd)));
maxScale = max(max(roiMapSum(:,pickROI,pickInd)));
scale = max([abs(minScale) abs(maxScale)]);
figure;
for ff=1:3
    subplot(2,3,ff);
    plotOnEgi(roiMapSum(:,pickROI,pickInd(ff)))
    caxis([-scale scale]);
    title(listROIs{pickROI})
end
pickROI = 18; pickInd = [46 45 4];
minScale = min(min(roiMapSum(:,pickROI,:)));
maxScale = max(max(roiMapSum(:,pickROI,:)));
scale = max([abs(minScale) abs(maxScale)]);
for ff=4:6
    subplot(2,3,ff);
    plotOnEgi(roiMapSum(:,pickROI,pickInd(ff-3)))
    caxis([-scale scale]);
    title(listROIs{pickROI})
end
set(gcf,'position',[100,100,1000,600])
saveas(gcf,'indROI/representative','png')

% for bilateral activation
figure;set(gcf,'position',[100,100,1200,800])
pickInd = [46 45 4];
pickROI = [1 2 1 2; 17 18 17 18; 1 2 17 18 ] ;
for tt = 1:size(pickROI,1)
    minScale = min(min(sum(roiMapSum(:,pickROI(tt,:),pickInd),2)));
    maxScale = max(max(sum(roiMapSum(:,pickROI(tt,:),pickInd),2)));
    scale = max([abs(minScale) abs(maxScale)]);
    for ff=1:3
        subplot(3,3,tt+3*(ff-1));
        plotOnEgi(sum(roiMapSum(:,pickROI(tt,:),pickInd(ff)),2))
        caxis([-scale scale]);
    end
end
saveas(gcf,'indROI/bilat3sbj','png')
figure;set(gcf,'position',[100,100,300,300])
data = mean(sum(roiMapSum(:,pickROI(tt,:),1:20),2),3);
plotOnEgi(data); caxis([-max(abs(data)) max(abs(data))]);
saveas(gcf,'indROI/bilatAverage','png')
    
figure;set(gcf,'position',[100,100,1000,800])
pickROI = [1 2; 15 16; 17 18] ;pickInd=1;
for tt = 1:size(pickROI,1)
    minScale = min(min(sum(roiMapSum(:,pickROI(tt,:),pickInd),2)));
    maxScale = max(max(sum(roiMapSum(:,pickROI(tt,:),pickInd),2)));
    scale = max([abs(minScale) abs(maxScale)]);
    subplot(3,1,tt);
    plotOnEgi(sum(roiMapSum(:,pickROI(tt,:),pickInd),2))
    caxis([-scale scale]);
end
saveas(gcf,'indROI/bilatS1','png')



% plot variability in amplitude
for rr=1:size(roiMapSum,2)
    for ss=1:size(roiMapSum,3)
        minAmp(rr,ss) = min(roiMapSum(:,rr,ss));  
        maxAmp(rr,ss) = max(roiMapSum(:,rr,ss));
    end
end
figure;
for rr=1:size(roiMapSum,2)
subplot(3,6,rr)
plot(minAmp(rr,:)); hold on; plot(maxAmp(rr,:))
title(listROIs(rr))
line(1:50,repmat(mean(minAmp(rr,:))- 3*std(minAmp(rr,:)),1,50),'Color','b')
line(1:50,repmat(mean(minAmp(rr,:))+ 3*std(minAmp(rr,:)),1,50),'Color','b')
line(1:50,repmat(mean(maxAmp(rr,:))- 3*std(maxAmp(rr,:)),1,50),'Color','r')
line(1:50,repmat(mean(maxAmp(rr,:))+ 3*std(maxAmp(rr,:)),1,50),'Color','r')
end


