function makeTemplate(projectName)
% creates a .mat with the avMap (elec x ROI) and list of ROIs
% e.g. run makeTemplate('Biosemi256')

%%% ATTENTION when plotting the topograhies: should be according to the system layout!!!

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated

listFiles = dir(['/Users/marleneponcet/Documents/data/skeriDATA/forward' projectName '/*.mat']);
if isempty(listFiles)
    disp('Cannot find the forward folder for this project. Please run makeFwd')
end
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
aa=0;
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
save(['averageMap' projectName '.mat'],'avMap','listROIs')

% plot if layout available
mm = round(max(max(abs(avMap))),-1);
projectName = lower(projectName);
if contains(projectName,'biosemi')
    figure('position', [200, 0, 1500, 800])
    colormap jet
    for roi=1:18
        subplot(3,6,roi)
        title(listROIs(roi))
        plotTopo(avMap(:,roi),[projectName '.lay'])
        caxis([-mm mm])
        % colorbar
    end
    saveas(gcf,['figures/averageMap' projectName],'png')
elseif strcmp(projectName,'egi256') ||  strcmp(projectName,'egi128') 
    figure('position', [200, 0, 1500, 800])
    for roi=1:18
        subplot(3,6,roi)
        title(listROIs(roi))
        plotOnEgi(avMap(:,roi)) % only for 256 & 128 electrodes
        caxis([-mm mm])
    end
    saveas(gcf,['figures/averageMap' projectName],'png')
end



