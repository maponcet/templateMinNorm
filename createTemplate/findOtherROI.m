%%% look in the forward models
% path to data
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128/';
dirList = dir([dataPath 'forward*']);
% list of ROIs in templates
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
% check which additional ROI is present 
for ff=1:length(dirList)
    load([dataPath dirList(ff).name])
    allROI = {roiInfo.name};
    index = contains(allROI,listROIs);
    % get ROIs that are not in the template list
    addROI{ff} = allROI(index == 0);
end

% get a count of how many sbj has each additional ROIs
clear countSbj newList
newList1 = addROI{20}; % sbj with max ROIs
newList2 = addROI{3};
newList3 = addROI{45};
newList = unique([newList1 newList2 newList3]);
for ff=1:length(addROI)
    countSbj(ff,:) = contains(newList,addROI{ff});
    if sum(countSbj(ff,:)) ~= length(addROI{ff})
        ff % check if sbj has additional ROI
    end
end
total = sum(countSbj);
otherROI = table(newList',total','VariableNames',{'ROI','nbOfSbj'});
writetable(otherROI)

% smallFOV NOT exact substraction of V1V-V1Vecc



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% look at non functional ROIs
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardEGI128_allROI/';
dirList = dir([dataPath 'forward*']);
otherROI = readtable('otherROI.txt');
listROIs = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};
functROI = [otherROI.ROI' listROIs];
% list of anatomical ROIs
for ff=1:length(dirList)
    load([dataPath dirList(ff).name])
    allROI = {roiInfo.name};
    index = contains(allROI,functROI);
    % get ROIs that are not in the template list
    addROI{ff} = allROI(index == 0);
end
% 2 sbj with insula in addition, all others have the same 70 anatomical
% regions
anatROI = addROI{1}(1:end-2); % remove the 2 unknown
% find overlapping indexes with visual ROIs
for ff=1:length(dirList)
    load([dataPath dirList(ff).name])
    % get indexes corresponding to visual areas
    index = find(contains({roiInfo.name},listROIs));
    allVisual = unique([roiInfo(index).meshIndices]);
    % count any overlap with anatomical ROIs
    indexAnat = find(contains({roiInfo.name},anatROI));
    for roiA = 1:length(indexAnat)
        totalMesh(ff,roiA) = length(roiInfo(indexAnat(roiA)).meshIndices);
        overlap(ff,roiA) = sum(ismember(roiInfo(indexAnat(roiA)).meshIndices,allVisual));
    end
end
% across sub
tt = sum(totalMesh);
oo = sum(overlap);
rr = sum(overlap) ./ sum(totalMesh);
mm = mean(totalMesh);

noOverlap = find(sum(overlap) == 0);
listExtROI = anatROI(noOverlap);
listExtROI'
% 1 ROI overlaps only in one hemisphere (3 indices!!). 'superiortemporal-L'
% Remove the other one 'superiortemporal-R'
listExtROI = listExtROI([1:end-7 end-5:end]);
nbMesh = mean(totalMesh(:,noOverlap([1:end-7 end-5:end])));
shortList = listExtROI(1:2:end);
nbMeshShort = nbMesh(1:2:end);
