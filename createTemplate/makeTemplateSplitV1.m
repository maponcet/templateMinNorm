function makeTemplateSplitV1(projectName)
% creates a structure templates with fields including electrode locations, 
% weights, list of ROIs 
% V1V and V1D instead of V1
% e.g. run makeTemplateSplitV1('EGI128')

% addpath if needed
if ~exist('ft_read_sens','file')
    addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork/fileio/
end


listFiles = dir(['/Users/marleneponcet/Documents/data/skeriDATA/forward' projectName '/*.mat']);
if isempty(listFiles)
    disp('Cannot find the forward folder for this project. Please run makeFwd')
end
listROIs = {'V1V-L', 'V1V-R', 'V1D-L', 'V1D-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
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

% % do a normalisation for each ROI -> unit norming
% % so that the betas represent microVolts (instead of microVolts/area size
% % as it is now)
% % unit norming is: all electrodes are squared and summed. These values are
% % then divided so that the total of the electrodes for each ROI (power) is
% % equal to 1
% regParam = sqrt(sum(avMap.^2,1));
% avMapNorm = bsxfun(@rdivide,avMap,regParam);

if contains(projectName,'standard_1005','IgnoreCase',true)
    templates = ft_read_sens('templates/standard_1005.elc');
elseif contains(projectName, 'EGI','IgnoreCase',true)
    templates = ft_read_sens(['templates/GSN-HydroCel-' projectName(4:end) '.sfp']);
end
templates.listROIs = listROIs;
templates.weights = avMap;
templates.reference = 0;
save(['templates' filesep 'template_V1Split_' projectName '.mat'],'templates')

