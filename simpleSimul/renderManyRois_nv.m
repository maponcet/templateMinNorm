
names = dir('/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/forward*');
for i=1:length(names)
    subjList{i} = names(i).name(end-12:end-4);
end

% subjList = {'skeri0001'}
% roiList = {'V1'}

roiList = {'V1', 'V2V', 'V2D', 'V3V','V3D', 'V4', 'V3A','LOC', 'MT'};

hemiName = {'-L' '-R'};
hemiLongName = {'left' 'right'};
cp(1,:)=[0.8   -0.27   -0.02];
cp(2,:)=[ -0.8   -0.27   -0.02];

ct=[0   -0.06    0.04];

visualView = [1:3];
infView = [4 6];
latView = [5 7:9];

for iSubj = 1:length(subjList)
    for iRoi = 1:length(roiList)
        
        roiDir = fullfile(getpref('mrCurrent','AnatomyFolder'),subjList{iSubj},'Standard','meshes','ROIs');
        roiLabel = ones(20484,1);
        
        for iHemi = 1:2
            roiFile = fullfile(roiDir,[roiList{iRoi} hemiName{iHemi} '.mat']);
            
            if ~exist(roiFile,'file')
                fprintf('file missing\n')
                continue;
            end
            load(roiFile);            
            roiLabel(ROI.meshIndices,:) = 2;
        end
        figure(100)
        clf;
        
        %         isRight=(1-(iHemi-1));
        %         thisHemiAlpha = [ isRight*ones(10242,1); (1-isRight)*ones(10242,1)];
        thisHemiAlpha = [ ones(10242,1); ones(10242,1)];
        
        cH = drawAnatomySurface(subjList{iSubj},'cortex');
        vertices = get(cH,'vertices');
        roiMean = mean(vertices(ROI.meshIndices,:),1);
        set(gca,'view',[0 0])
        %         if ismember(iRoi,infView)
        %             set(gca,'view',[180 270])
        %         elseif ismember(iRoi,latView)
        %             set(gca,'view',[(2*iHemi-3)*90 0])
        %         elseif ismember(iRoi,visualView)
        %             cp(iHemi,2:3) = roiMean(2:3);
        %             campos(cp(iHemi,:));
        %             camtarget(roiMean);%camtarget(ct);
        %         end
        set(cH,'faceVertexAlphaData',thisHemiAlpha,'facealpha',0.75, ...
            'facevertexCData',roiLabel,'CDataMapping','direct','facecolor','interp');
        colormap( [.8 .8 .8;1 0 0]);
        material([.3 .5 0 .2]);
        %    light('position',[0 .707 .707],'style','infinite')
        %    light('position',[0 -.707 .707],'style','infinite')
        light('position',[0 0 1],'style','local')
        lighting phong
        set(gcf,'color','w')
        camproj('orthographic');grid minor;
        th=text(roiMean(1),roiMean(2),roiMean(3)-.03,[subjList{iSubj} ' ' hemiLongName{iHemi} ' hemisphere']);
        set(th,'fontsize',20);
        set(gcf,'position',[54 179 900 747]);
        ylabel('Distance in Meters')
        filename = ['renderROI/' 'withGrid_' subjList{iSubj} '_' roiList{iRoi} '.png'];
        saveas(gcf,filename);
        
    end
end

