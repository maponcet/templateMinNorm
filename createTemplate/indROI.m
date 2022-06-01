% plot different representative response
% using 256 electrodes and other plotting functions
% Use forward matrices from specified folder (can be any montage as long as
% the electrode layout matches)
% Addition: do for 128 EGI

clearvars
% add path for plotTopo  ft_prepare_layout plotOnEGI
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/svndl_code/alesToolbox


%%% 256 EGI
listFiles = dir('/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI256/*.mat');
elecLayout = 'layout/GSN-HydroCel-256.sfp'; % 'layout/biosemi256.lay'

%%% 128 EGI
listFiles = dir('/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/*.mat');
elecLayout = 'layout/GSN-HydroCel-128.sfp'; % 'layout/biosemi256.lay'

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

pickROI = 2; pickInd = [2 30 19];
minScale = min(min(roiMapSum(:,pickROI,pickInd)));
maxScale = max(max(roiMapSum(:,pickROI,pickInd)));
scale = max([abs(minScale) abs(maxScale)]);
figure;
for ff=1:3
    subplot(2,3,ff);
    plotTopo(roiMapSum(:,pickROI,pickInd(ff)),elecLayout)
    caxis([-scale scale]);
    title(listROIs{pickROI})
    colorcet('D1')
end
pickROI = 18; pickInd = [46 45 4];
minScale = min(min(roiMapSum(:,pickROI,:)));
maxScale = max(max(roiMapSum(:,pickROI,:)));
scale = max([abs(minScale) abs(maxScale)]);
for ff=4:6
    subplot(2,3,ff);
    plotTopo(roiMapSum(:,pickROI,pickInd(ff-3)),elecLayout)
    caxis([-scale scale]);
    title(listROIs{pickROI})
    colorcet('D1')
end
set(gcf,'position',[100,100,1000,600])
saveas(gcf,'indROI/representative128','png')
saveas(gcf,'indROI/representative128','fig')
print(gcf,'indROI/representative128','-depsc')

% same but EEGlab topoplot
pickROI = 2; pickInd = [2 30 19];
minScale = min(min(roiMapSum(:,pickROI,pickInd)));
maxScale = max(max(roiMapSum(:,pickROI,pickInd)));
scale = max([abs(minScale) abs(maxScale)]);
figure;
for ff=1:3
    subplot(2,3,ff);
    topoplot(roiMapSum(:,pickROI,pickInd(ff)),elecLayout,'colormap',colorcet('D1'),'electrodes','on' )
    caxis([-scale scale]);
    title(listROIs{pickROI})
    colorcet('D1')
end
pickROI = 18; pickInd = [46 45 4];
minScale = min(min(roiMapSum(:,pickROI,:)));
maxScale = max(max(roiMapSum(:,pickROI,:)));
scale = max([abs(minScale) abs(maxScale)]);
for ff=4:6
    subplot(2,3,ff);
    topoplot(roiMapSum(:,pickROI,pickInd(ff-3)),elecLayout,'colormap',colorcet('D1'),'electrodes','on' )
    caxis([-scale scale]);
    title(listROIs{pickROI})
    colorcet('D1')
end
set(gcf,'position',[100,100,1000,600])
saveas(gcf,'indROI/representative256topoplot','png')



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
        plotTopo(sum(roiMapSum(:,pickROI(tt,:),pickInd(ff)),2),elecLayout)
        caxis([-scale scale]);
        colorcet('D1')
    end
end
saveas(gcf,'indROI/bilat3sbj_256','png')

figure;set(gcf,'position',[100,100,300,300])
data = mean(sum(roiMapSum(:,pickROI(tt,:),1:20),2),3);
plotTopo(data,elecLayout); caxis([-max(abs(data)) max(abs(data))]);colorcet('D1')
saveas(gcf,'indROI/bilatAverage_256','png')
    
figure;set(gcf,'position',[100,100,1000,800])
pickROI = [1 2; 15 16; 17 18] ;pickInd=1;
for tt = 1:size(pickROI,1)
    minScale = min(min(sum(roiMapSum(:,pickROI(tt,:),pickInd),2)));
    maxScale = max(max(sum(roiMapSum(:,pickROI(tt,:),pickInd),2)));
    scale = max([abs(minScale) abs(maxScale)]);
    subplot(3,1,tt);
    plotTopo(sum(roiMapSum(:,pickROI(tt,:),pickInd),2),elecLayout)
    caxis([-scale scale]);colorcet('D1');
end
saveas(gcf,'indROI/bilatS1_256','png')






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
    colorcet('D1')
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
    colorcet('D1')
end
set(gcf,'position',[100,100,1000,600])
saveas(gcf,'indROI/representative128plotOnEgi','png')
saveas(gcf,'indROI/representative128plotOnEgi','fig')
print(gcf,'indROI/representative128plotOnEgi','-depsc')
% same but EEGlab topoplot
pickROI = 2; pickInd = [2 30 19];
minScale = min(min(roiMapSum(:,pickROI,pickInd)));
maxScale = max(max(roiMapSum(:,pickROI,pickInd)));
scale = max([abs(minScale) abs(maxScale)]);
figure;
for ff=1:3
    subplot(2,3,ff);
    topoplot(roiMapSum(:,pickROI,pickInd(ff)),elecLayout,'colormap',colorcet('D1'),'electrodes','on' )
    caxis([-scale scale]);
    title(listROIs{pickROI})
    colorcet('D1')
end
pickROI = 18; pickInd = [46 45 4];
minScale = min(min(roiMapSum(:,pickROI,:)));
maxScale = max(max(roiMapSum(:,pickROI,:)));
scale = max([abs(minScale) abs(maxScale)]);
for ff=4:6
    subplot(2,3,ff);
    topoplot(roiMapSum(:,pickROI,pickInd(ff-3)),elecLayout,'colormap',colorcet('D1'),'electrodes','on' )
    caxis([-scale scale]);
    title(listROIs{pickROI})
    colorcet('D1')
end
set(gcf,'position',[100,100,1000,600])
saveas(gcf,'indROI/representative128eeglab','png')
saveas(gcf,'indROI/representative128eeglab','fig')
print(gcf,'indROI/representative128eeglab','-depsc')


pickROI = 2; pickInd = [2 30 19];
minScale = min(min(roiMapSum(:,pickROI,pickInd)));
maxScale = max(max(roiMapSum(:,pickROI,pickInd)));
scale = max([abs(minScale) abs(maxScale)]);
figure;
for ff=1:3
    subplot(2,3,ff);
    plotContourOnScalp(roiMapSum(:,pickROI,pickInd(ff)),'skeri0044',...
        '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
    view(20,35);camproj('perspective');axis off
    caxis([-scale scale]);title(listROIs{pickROI})
    title(listROIs{pickROI})
    colorcet('D1')
end
pickROI = 18; pickInd = [46 45 4];
minScale = min(min(roiMapSum(:,pickROI,:)));
maxScale = max(max(roiMapSum(:,pickROI,:)));
scale = max([abs(minScale) abs(maxScale)]);
for ff=4:6
    subplot(2,3,ff);
    plotContourOnScalp(roiMapSum(:,pickROI,pickInd(ff-3)),'skeri0044',...
        '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
    view(20,35);camproj('perspective');axis off
    caxis([-scale scale]);title(listROIs{pickROI})
    colorcet('D1')
end
set(gcf,'position',[100,100,800,800])
saveas(gcf,'indROI/representative3D','png')
saveas(gcf,'indROI/representative3D','fig')
print('indROI/representative3D','-depsc')
