function makeTemplate(projectName)
% creates a .mat with the avMap (elec x ROI) and list of ROIs
% e.g. run makeTemplate('Biosemi256')

%%% ATTENTION when plotting the topograhies: should be according to the system layout!!!

addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/svndl_code/alesToolbox

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

% do a normalisation for each ROI -> unit norming
% so that the betas represent microVolts (instead of microVolts/area size
% as it is now)
% unit norming is: all electrodes are squared and summed. These values are
% then divided so that the total of the electrodes for each ROI (power) is
% equal to 1
regParam = sqrt(sum(avMap.^2,1));
avMapNorm = bsxfun(@rdivide,avMap,regParam);

if contains(projectName,'standard_1005','IgnoreCase',true)
    templates = ft_read_sens('templates/standard_1005.elc');
elseif contains(projectName, 'EGI','IgnoreCase',true)
    templates = ft_read_sens(['templates/GSN-HydroCel-' projectName(4:end) '.sfp']);
end
templates.listROIs = listROIs;
templates.weights = avMap;
templates.reference = 0;
save(['templates' filesep 'template_' projectName '.mat'],'templates')

% plot if layout available
mm = round(max(max(abs(avMap))),-1);
% nn = round(max(max(abs(avMapNorm))),1);
loc = [1:9;10:18]; loc = loc(:);
if contains(projectName,'biosemi') 
    figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        title(listROIs(roi))
        plotTopo(avMap(:,roi),['layout/' projectName '.lay'])
        caxis([-mm mm])
        colorcet('D1') 
    end
    saveas(gcf,['figures/averageMap' projectName 'plotTopo'],'png')
%     figure('position', [200, 1000, 2000, 500])
%     for roi=1:18
%         subplot(2,9,loc(roi))
%         title(listROIs(roi))
%         plotTopo(avMapNorm(:,roi),['layout/' projectName '.lay'])
%         caxis([-nn nn])
%         colorcet('D1') 
%     end
%     saveas(gcf,['figures/averageMapNorm' projectName],'png')
    
%%%%%%
elseif contains(projectName,'EGI')
    nbChan = sscanf(projectName,'EGI%d');
    figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        plotTopo(avMap(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp']);
        caxis([-mm mm]);title(listROIs(roi));
        colorcet('D1') 
    end
    saveas(gcf,['figures/averageMap' projectName 'plotTopo'],'png')
    figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        topoplot(avMap(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp'],'colormap',colorcet('D1'),'electrodes','on' );
        caxis([-mm mm]);title(listROIs(roi));
    end
    saveas(gcf,['figures/averageMap' projectName 'EEGlab'],'png')
    
    if nbChan > 64 
        figure('position', [200, 1000, 2000, 500])
        for roi=1:18
            subplot(2,9,loc(roi))
            title(listROIs(roi))
            plotOnEgi(avMap(:,roi)) % only for 256 & 128 electrodes
            caxis([-mm mm])
            colorcet('D1') 
        end
        saveas(gcf,['figures/averageMap' projectName 'plotOnEgi'],'png')
        
%         figure('position', [200, 1000, 2000, 500])
%         for roi=1:18
%             subplot(2,9,loc(roi))
%             title(listROIs(roi))
%             plotOnEgi(avMapNorm(:,roi)) % only for 256 & 128 electrodes
%             caxis([-nn nn])
%             colorcet('D1') 
%         end
%         saveas(gcf,['figures/averageMapNorm' projectName],'png')
    end
    
%%%%%%%%%     
elseif strcmp(projectName,'standard_1005')
    figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        plotTopo(avMap(:,roi),'layout/Standard-10-5-Cap385.sfp');
        caxis([-mm mm]);title(listROIs(roi));
        colorcet('D1') 
    end
    saveas(gcf,['figures/averageMap' projectName],'png')
    saveas(gcf,['figures/averageMap' projectName],'fig')
    print(gcf,['figures/averageMap' projectName],'-depsc')
    figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        % although same file for the layout, topoplot use readlocs from
        % EEGlab and converts EGI Cartesian coordinates to Matlab/EEGLAB xyz coordinates.
        topoplot(avMap(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','on' );
        caxis([-mm mm]);title(listROIs(roi))
    end
    saveas(gcf,['figures/averageMap' projectName 'EEGlab'],'png')
    figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        topoplot(avMap(:,roi),'Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','off' );
        caxis([-mm mm]);title(listROIs(roi));
    end
    saveas(gcf,['figures/averageMap' projectName 'EEGlabNoElec'],'png')
end



% % fieldtrip layout gets error.. :/
% dirlist  = dir('fieldtrip/*.*');
% filename = {dirlist(~[dirlist.isdir]).name}';
% cfg = [];
% cfg.layout = 'fieldtrip/biosemi64.lay';
% layout = ft_prepare_layout(cfg);
% figure; ft_plot_topo(cfg,avMap(:,roi));
% figure; ft_topoplotER(cfg,abs(avMap(:,roi)));
% ft_layoutplot(layout);

%  cfg.layout = 'EEG1005_modif.lay';
%  ft_prepare_layout(cfg);
%  figure
%  ft_layoutplot(cfg);
% figure; plotTopo(avMap(4:end-7,1),'EEG1005_modif.lay');
