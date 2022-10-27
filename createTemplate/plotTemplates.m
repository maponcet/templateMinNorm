%%%% plot templates for the different montages
% could use colorcet('D9') for a more white 0, lighter colours
addpath('/Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated');
addpath(genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork'));
addpath(genpath('/Users/marleneponcet/Documents/Git/templateMinNorm/simpleSimul/subfunctions'));
addpath('/Users/marleneponcet/Documents/Git/svndl_code/alesToolbox');

listTemplates = dir('templates/template*');

%%% 3D only right hemisphere (works only for EGI128)
load('templates/template_EGI128.mat')
listROIs = templates.listROIs;
numROI = length(listROIs);
%%% left
figure; pict = 1;
for roi=1:2:numROI 
    subplot(1,numROI/2,pict)
    title(listROIs(roi))
    plotContourOnScalp(templates.weights(:,roi),'skeri0044',...
        '/Users/marleneponcet/Documents/Git/templateMinNorm/simpleSimul/realData/eegdata/')
    view(-30,35)
    camproj('perspective')
    axis off
    pict = pict+1;
end
set(gcf,'position', [200, 0, 2000, 300])
saveas(gcf,'figures/averageMap50_3DturnL','png')
saveas(gcf,'figures/averageMap50_3DturnL','fig')
print(gcf,'figures/averageMap50_3DturnL','-depsc')

figure; pict = 1;
for roi=2:2:numROI 
    subplot(1,numROI/2,pict)
    title(listROIs(roi))
    plotContourOnScalp(templates.weights(:,roi),'skeri0044',...
        '/Users/marleneponcet/Documents/Git/templateMinNorm/simpleSimul/realData/eegdata/')
    view(30,35)
    camproj('perspective')
    axis off
    pict = pict+1;
end
set(gcf,'position', [200, 0, 2000, 300])
saveas(gcf,'figures/averageMap50_3DturnR','png')
saveas(gcf,'figures/averageMap50_3DturnR','fig')
print(gcf,'figures/averageMap50_3DturnR','-depsc')

%%% 2D
% location ROI for subplot

for tt=1:length(listTemplates)
    clear avMap
    load([listTemplates(tt).folder '/' listTemplates(tt).name])
    listROIs = templates.listROIs;
    numROI = length(listROIs);
    loc = [1:numROI/2;numROI/2+1:numROI]; loc = loc(:);
    ttName = listTemplates(tt).name(10:end-4);
    mm = round(max(max(abs(templates.weights))),-1);
    
    if contains(ttName,'Biosemi')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:numROI
            subplot(2,numROI/2,loc(roi))
            title(listROIs(roi))
            plotTopo(templates.weights(:,roi),['layout/' ttName '.lay'])
            caxis([-mm mm])
            colorcet('D1')
        end
        saveas(gcf,['figures/averageMap' ttName 'plotTopo'],'png')
        saveas(gcf,['figures/averageMap' ttName 'plotTopo'],'fig')
        print(['figures/averageMap' ttName 'plotTopo'],'-depsc')

        %%%%%%
    elseif contains(ttName,'EGI')
%         nbChan = sscanf(ttName,'%*sEGI%d');
        nbChan = extractAfter(ttName,'EGI');
        figure('position', [200, 1000, 2000, 500])
        for roi=1:numROI
            subplot(2,numROI/2,loc(roi))
            plotTopo(templates.weights(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp']);
            caxis([-mm mm]);title(listROIs(roi));
            colorcet('D1')
        end
        saveas(gcf,['figures/averageMap' ttName 'plotTopo'],'png')
        saveas(gcf,['figures/averageMap' ttName 'plotTopo'],'fig')
        print(gcf,['figures/averageMap' ttName 'plotTopo'],'-depsc')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:numROI
            subplot(2,numROI/2,loc(roi))
            topoplot(templates.weights(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp'],'colormap',colorcet('D1'),'electrodes','on' );
            caxis([-mm mm]);title(listROIs(roi));
        end
        saveas(gcf,['figures/averageMap' ttName 'EEGlab'],'png')
        saveas(gcf,['figures/averageMap' ttName 'EEGlab'],'fig')
        print(gcf,['figures/averageMap' ttName 'EEGlab'],'-depsc')
        
        %%%%% plotOnEgi only for 256 & 128 electrodes
        if nbChan > 64
            figure('position', [200, 1000, 2000, 500])
            for roi=1:numROI
                subplot(2,numROI/2,loc(roi))
                title(listROIs(roi))
                plotOnEgi(templates.weights(:,roi)) 
                caxis([-mm mm])
                colorcet('D1')
            end
            saveas(gcf,['figures/averageMap' ttName 'plotOnEgi'],'png')
            saveas(gcf,['figures/averageMap' ttName 'plotOnEgi'],'fig')
            print(gcf,['figures/averageMap' ttName 'plotOnEgi'],'-depsc')
        end
        
        %%%%%%%%%
    elseif strcmp(ttName,'Standard_1005')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:numROI
            subplot(2,numROI/2,loc(roi))
            plotTopo(templates.weights(:,roi),'layout/Standard-10-5-Cap385.sfp');
            caxis([-mm mm]);title(listROIs(roi));
            colorcet('D1')
        end
        saveas(gcf,['figures/averageMap' ttName],'png')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:numROI
            subplot(2,numROI/2,loc(roi))
            topoplot(templates.weights(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','on' );
            caxis([-mm mm]);title(listROIs(roi));
        end
        saveas(gcf,['figures/averageMap' ttName 'EEGlab'],'png')
        saveas(gcf,['figures/averageMap' ttName 'EEGlab'],'fig')
        print(gcf,['figures/averageMap' ttName 'EEGlab'],'-depsc')
        figure('position', [200, 1000, 2000, 500])
        for roi=1:numROI
            subplot(2,numROI/2,loc(roi))
            topoplot(templates.weights(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','off' );
            caxis([-mm mm]);title(listROIs(roi));
        end
        saveas(gcf,['figures/averageMap' ttName 'EEGlabNoElec'],'png')
    end
    
end


% %%%% plot layout
% figure; topoplot([],['layout/GSN-HydroCel-' num2str(nbChan) '.sfp'],'colormap',colorcet('D1'),'electrodes','labels' );
% figure; topoplot([],['layout/Standard-10-5-Cap385.sfp'],'colormap',colorcet('D1'),'electrodes','labels' );
% 
% listLayout = dir('layout/');
% listLayout = listLayout([listLayout.bytes]>0); % remove files with no data (.)
% for tt=1:length(listLayout)
%     ttName = listLayout(tt).name;
%     cfg=[];
%     cfg.layout = ['layout/' ttName];
%     layout = ft_prepare_layout(cfg);
%     cfg = [];
%     cfg.layout = layout;
%     ft_layoutplot(cfg);title(ttName)
%     saveas(gcf,['layout/layout' ttName(1:end-4)],'fig')
% end
    