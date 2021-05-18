% count the nb of overlapping 

dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
nbSbj = length(dirList);

listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};

for iSubj=1:nbSbj
    clear fwdMatrix roiInfo idxROI
    % fwd file
    load([dataPath dirList(iSubj).name])
    % indexes of the common ROIs  
    idxROI = find(ismember({roiInfo.name},listROIs));
    % check that all rois are there
    if length(idxROI) ~= length(listROIs)
        fprintf('nb of ROIs mismatch for %s \n',dirList(iSubj).name)
    else
        % go through each ROI and compare mesh indices 
        for rr=1:length(idxROI)
            curROI = find(strcmp(listROIs(rr),{roiInfo.name}));
            otherIdx = unique([roiInfo(setdiff(idxROI,curROI)).meshIndices]);
            currIdx = roiInfo(curROI).meshIndices;
            overlap(iSubj,rr) = sum(ismember(currIdx,otherIdx));
            total(iSubj,rr) = length(currIdx);
            prctOverlap(iSubj,rr) = overlap(iSubj,rr)/total(iSubj,rr)*100;
        end
    end
end

figure;hold on;
bar(mean(prctOverlap))
ylabel('% overlap with any other 17 ROIs')
xticks([1:18])
xticklabels(listROIs)
saveas(gcf,['figures' filesep 'meshOverlap'],'png')
