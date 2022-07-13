
%%% topography 250 ms (as in Lim)
load('realDataOutput/realData.mat')
load('realDataOutput/realDataLasso.mat')
load('realDataOutput/timeline.mat')

timeToPlot = find(round(t)==250); % =250ms
addpath('/Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated/')
if ~exist('ft_prepare_layout','file')
    addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork/
end
if ~exist('plotContourOnScalp','file')
    addpath(genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork'));
end


layoutData = 'realData/GSN-HydroCel-128.sfp';
mm = round(max([abs(mean(Ypca(:,:,timeToPlot),1))'; abs(reTemplates(:,timeToPlot)); ...
    abs(reWhole(:,timeToPlot)); abs(reROI(:,timeToPlot)); abs(reLasso(:,timeToPlot))]),1);
% 2D fieldtrip
figure;
subplot(2,3,1);plotTopo(mean(Ypca(:,:,timeToPlot),1),layoutData);caxis([-mm mm]);title('data')
subplot(2,3,2);plotTopo(reTemplates(:,timeToPlot),layoutData);caxis([-mm mm]);title('Template')
subplot(2,3,3);plotTopo(reROI(:,timeToPlot),layoutData);caxis([-mm mm]);title('Individual')
subplot(2,3,4);plotTopo(reWhole(:,timeToPlot),layoutData);caxis([-mm mm]);title('Individual Subset')
subplot(2,3,5);plotTopo(reLasso(:,timeToPlot),layoutData);caxis([-mm mm]);title('Lasso')
colorcet('D1')
set(gcf, 'Position', [100 500 1000 600]);
saveas(gcf,['figures/realDataTopo2D'],'png')
saveas(gcf,['figures/realDataTopo2D'],'fig')
print(gcf,['figures/realDataTopo2D'],'-depsc')


% dd = max(abs(mean(Ypca(:,:,timeToPlot),1)));
% tt = max(abs(reTemplates(:,timeToPlot)));
% ww = max(abs(reWhole(:,timeToPlot)));
% rr = max(abs(reROI(:,timeToPlot)));
% ll = max(abs(reLasso(:,timeToPlot)));
% figure;
% subplot(2,3,1);plotTopo(mean(Ypca(:,:,timeToPlot),1),layoutData);caxis([-dd dd]);title('data')
% subplot(2,3,2);plotTopo(reTemplates(:,timeToPlot),layoutData);caxis([-tt tt]);title('Template')
% subplot(2,3,3);plotTopo(reWhole(:,timeToPlot),layoutData);caxis([-ww ww]);title('Individual')
% subplot(2,3,4);plotTopo(reROI(:,timeToPlot),layoutData);caxis([-rr rr]);title('Individual Subset')
% subplot(2,3,5);plotTopo(reLasso(:,timeToPlot),layoutData);caxis([-ll ll]);title('Lasso')
% colorcet('D1')
% set(gcf, 'Position', [100 500 1000 600]);
% saveas(gcf,['figures/realDataTopo2D'],'png')
% saveas(gcf,['figures/realDataTopo2D'],'fig')
% print(gcf,['figures/realDataTopo2D'],'-depsc')


% 3D
figure;
subplot(2,3,1);
plotContourOnScalp(mean(Ypca(:,:,timeToPlot),1),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('data')
subplot(2,3,2);
plotContourOnScalp(reTemplates(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('template')
subplot(2,3,3);
plotContourOnScalp(reROI(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('constrained')
subplot(2,3,4);
plotContourOnScalp(reWhole(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('tailored')
subplot(2,3,5);
plotContourOnScalp(reLasso(:,timeToPlot),'skeri0044',...
    '/Users/marleneponcet/Documents/Git/templateMinNorm/PlosOne/github-archive/datafiles/eegdata/')
set(gcf, 'Position', [100 500 1000 600]);
view(20,35);camproj('perspective');axis off;caxis([-mm mm]);title('lasso')
saveas(gcf,['figures/realDataTopo3D'],'png')
saveas(gcf,['figures/realDataTopo3D'],'fig')
print(gcf,['figures/realDataTopo3D'],'-depsc')