%%%% plot templates for the different montages
% could use colorcet('D9') for a more white 0, lighter colours

listTemplates = dir('templates/template*');

%%% 3D only right hemisphere (works only for EGI128)
load([listTemplates(5).folder '/' listTemplates(5).name])
addpath(genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork'));
figure; pict = 1;
for roi=2:2:18
    subplot(1,9,pict)
    title(listROIs(roi))
    plotContourOnScalp(avMap(:,roi),'skeri0044',...
        '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
    pict = pict+1;
end
set(gcf,'position', [200, 0, 2000, 300])
saveas(gcf,'averageMap50_3Dturn','png')


%%% 2D
% location ROI for subplot
loc = [1:9;10:18]; loc = loc(:);

for tt=1:length(listTemplates)
    clear avMap
    load([listTemplates(tt).folder '/' listTemplates(tt).name])
    ttName = listTemplates(tt).name(10:end-4);
    mm = round(max(max(abs(avMap))),-1);
    
    if contains(ttName,'Biosemi')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:18
            subplot(2,9,loc(roi))
            title(listROIs(roi))
            plotTopo(avMap(:,roi),['layout/' ttName '.lay'])
            caxis([-mm mm])
            colorcet('D1')
        end
        saveas(gcf,['figures/averageMap' ttName 'plotTopo'],'png')
        
        %%%%%%
    elseif contains(ttName,'EGI')
        nbChan = sscanf(ttName,'EGI%d');
        figure('position', [200, 1000, 2000, 500])
        for roi=1:18
            subplot(2,9,loc(roi))
            plotTopo(avMap(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp']);
            caxis([-mm mm]);title(listROIs(roi));
            colorcet('D1')
        end
        saveas(gcf,['figures/averageMap' ttName 'plotTopo'],'png')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:18
            subplot(2,9,loc(roi))
            topoplot(avMap(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp'],'colormap',colorcet('D1'),'electrodes','on' );
            caxis([-mm mm]);title(listROIs(roi));
        end
        saveas(gcf,['figures/averageMap' ttName 'EEGlab'],'png')
        
%         %%%%% plotOnEgi: maps are flipped for 256 electrodes!! 
%         if nbChan > 64
%             figure('position', [200, 1000, 2000, 500])
%             for roi=1:18
%                 subplot(2,9,loc(roi))
%                 title(listROIs(roi))
%                 plotOnEgi(avMap(:,roi)) % only for 256 & 128 electrodes
%                 caxis([-mm mm])
%                 colorcet('D1')
%             end
%             saveas(gcf,['figures/averageMap' ttName 'plotOnEgi'],'png')
%         end
        
        %%%%%%%%%
    elseif strcmp(ttName,'Standard_1005')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:18
            subplot(2,9,loc(roi))
            plotTopo(avMap(:,roi),'layout/Standard-10-5-Cap385.sfp');
            caxis([-mm mm]);title(listROIs(roi));
            colorcet('D1')
        end
        saveas(gcf,['figures/averageMap' ttName],'png')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:18
            subplot(2,9,loc(roi))
            topoplot(avMap(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','on' );
            caxis([-mm mm]);title(listROIs(roi));
        end
        saveas(gcf,['figures/averageMap' ttName 'EEGlab'],'png')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:18
            subplot(2,9,loc(roi))
            topoplot(avMap(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','off' );
            caxis([-mm mm]);title(listROIs(roi));
        end
        saveas(gcf,['figures/averageMap' ttName 'EEGlabNoElec'],'png')
    end
    
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

%%%% plot layout
figure; topoplot([],['layout/GSN-HydroCel-' num2str(nbChan) '.sfp'],'colormap',colorcet('D1'),'electrodes','labels' );
figure; topoplot([],['layout/Standard-10-5-Cap385.sfp'],'colormap',colorcet('D1'),'electrodes','labels' );

% for the plotTopo 
listLayout = dir('layout/');
listLayout = listLayout([listLayout.bytes]>0); % remove files with no data (.)
for tt=1:length(listLayout)
    ttName = listLayout(tt).name;
    cfg=[];
    cfg.layout = ['layout/' ttName];
    layout = ft_prepare_layout(cfg);
    cfg = [];
    cfg.layout = layout;
    figure;ft_layoutplot(cfg);
    saveas(gcf,['figures/layout' ttName(1:end-4)],'png')
end
    