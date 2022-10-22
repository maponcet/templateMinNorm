%%%% plot templates for the different montages
% could use colorcet('D9') for a more white 0, lighter colours
addpath('/Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated');
addpath(genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork'));
addpath(genpath('/Users/marleneponcet/Documents/Git/templateMinNorm/simpleSimul/subfunctions'));
addpath('/Users/marleneponcet/Documents/Git/svndl_code/alesToolbox');

listTemplates = dir('templates/template*');


%%% 3D only right hemisphere (works only for EGI128)
load('templates/template_EGI128.mat')
roiNames = {'V1','V2','V3','V4','V3A','LOC','MT'};
% which ROI to pool (bilateral + ventral/dorsal)
poolROI{1} = [1 2];
poolROI{2} = 3:6;
poolROI{3} = 7:10;
poolROI{4} = 11:12;
poolROI{5} = 13:14;
poolROI{6} = 15:16;
poolROI{7} = 17:18;

% plot
figure; pict = 1;
for roi=1:length(roiNames)
    subplot(1,length(roiNames),roi)
    title(roiNames(roi))
    plotContourOnScalp(sum(templates.weights(:,poolROI{roi}),2),'skeri0044',...
        '/Users/marleneponcet/Documents/Git/templateMinNorm/simpleSimul/realData/eegdata/')
%     view(-30,35)
    camproj('perspective')
    axis off
end
set(gcf,'position', [200, 0, 2000, 300])
saveas(gcf,'figures/fullField_3D','png')
saveas(gcf,'figures/fullField_3D','fig')
print(gcf,'figures/fullField_3D','-depsc')

%%% 2D
for tt=1:length(listTemplates)
    clear avMap bilatWeight
    load([listTemplates(tt).folder '/' listTemplates(tt).name])
    ttName = listTemplates(tt).name(10:end-4);
    
    % do the sum before plotting to pick the right scale
    for roi=1:length(roiNames)
        bilatWeight(:,roi) = sum(templates.weights(:,poolROI{roi}),2);
    end
    mm = round(max(max(abs(bilatWeight))),-1);
    
    if contains(ttName,'Biosemi')
        figure('position', [200, 1000, 2000, 250])
        for roi=1:length(roiNames)
            subplot(1,length(roiNames),roi)
            title(roiNames(roi))
            plotTopo(bilatWeight(:,roi),['layout/' ttName '.lay'])
            caxis([-mm mm])
            colorcet('D1')
        end
        saveas(gcf,['figures/fullField' ttName 'plotTopo'],'png')
        saveas(gcf,['figures/fullField' ttName 'plotTopo'],'fig')
        print(['figures/fullField' ttName 'plotTopo'],'-depsc')

        %%%%%%
    elseif contains(ttName,'EGI')
        nbChan = sscanf(ttName,'EGI%d');
        figure('position', [200, 1000, 2000, 250])
        for roi=1:length(roiNames)
            subplot(1,length(roiNames),roi)
            plotTopo(bilatWeight(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp']);
            caxis([-mm mm]);title(roiNames(roi));
            colorcet('D1')
        end
        saveas(gcf,['figures/fullField' ttName 'plotTopo'],'png')
        saveas(gcf,['figures/fullField' ttName 'plotTopo'],'fig')
        print(gcf,['figures/fullField' ttName 'plotTopo'],'-depsc')
        figure('position', [200, 1000, 2000, 250])
        for roi=1:length(roiNames)
            subplot(1,length(roiNames),roi)
            topoplot(bilatWeight(:,roi),['layout/GSN-HydroCel-' num2str(nbChan) '.sfp'],'colormap',colorcet('D1'),'electrodes','on' );
            caxis([-mm mm]);title(roiNames(roi));
        end
        saveas(gcf,['figures/fullField' ttName 'EEGlab'],'png')
        saveas(gcf,['figures/fullField' ttName 'EEGlab'],'fig')
        print(gcf,['figures/fullField' ttName 'EEGlab'],'-depsc')
        
        %%%%% plotOnEgi only for 256 & 128 electrodes
        if nbChan > 64
            figure('position', [200, 1000, 2000, 250])
            for roi=1:length(roiNames)
                subplot(1,length(roiNames),roi)
                title(roiNames(roi))
                plotOnEgi(bilatWeight(:,roi)) 
                caxis([-mm mm])
                colorcet('D1')
            end
            saveas(gcf,['figures/fullField' ttName 'plotOnEgi'],'png')
            saveas(gcf,['figures/fullField' ttName 'plotOnEgi'],'fig')
            print(gcf,['figures/fullField' ttName 'plotOnEgi'],'-depsc')
        end
        
        %%%%%%%%%
    elseif strcmp(ttName,'Standard_1005')
        figure('position', [200, 1000, 2000, 250])
        for roi=1:length(roiNames)
            subplot(1,length(roiNames),roi)
            plotTopo(bilatWeight(:,roi),'layout/Standard-10-5-Cap385.sfp');
            caxis([-mm mm]);title(roiNames(roi));
            colorcet('D1')
        end
        saveas(gcf,['figures/fullField' ttName],'png')
        figure('position', [200, 1000, 2000, 250])
        for roi=1:length(roiNames)
            subplot(1,length(roiNames),roi)
            topoplot(bilatWeight(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','on' );
            caxis([-mm mm]);title(roiNames(roi));
        end
        saveas(gcf,['figures/fullField' ttName 'EEGlab'],'png')
        saveas(gcf,['figures/fullField' ttName 'EEGlab'],'fig')
        print(gcf,['figures/fullField' ttName 'EEGlab'],'-depsc')
        figure('position', [200, 1000, 2000, 250])
        for roi=1:length(roiNames)
            subplot(1,length(roiNames),roi)
            topoplot(bilatWeight(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','off' );
            caxis([-mm mm]);title(roiNames(roi));
        end
        saveas(gcf,['figures/fullField' ttName 'EEGlabNoElec'],'png')
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
    