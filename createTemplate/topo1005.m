addpath /Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork
addpath /Users/marleneponcet/Documents/Git/ssvepTesting/biosemiUpdated
addpath /Users/marleneponcet/Documents/Git/svndl_code/alesToolbox

load('templates/template_Standard_1005.mat')

mm = round(max(max(abs(avMap))),-1);
loc = [1:9;10:18]; loc = loc(:);
figure('position', [200, 1000, 2000, 500])
for roi=1:18
    subplot(2,9,loc(roi))
    plotTopo(avMap(:,roi),'layout/Standard-10-5-Cap385.sfp');
    caxis([-mm mm]);title(listROIs(roi));
    colorcet('D1')
end

% bilateral activation
bilat = cell2mat(arrayfun(@(x) sum(avMap(:,[x x+1]),2),1:2:18,'UniformOutput',0));
listBilat = listROIs(1:2:18);
mm = round(max(max(abs(bilat))),-1);
figure('position', [200, 1000, 2000, 250])
for roi=1:9
    subplot(1,9,roi)
    plotTopo(bilat(:,roi),'layout/Standard-10-5-Cap385.sfp');
    caxis([-mm mm]);title(listBilat(roi));
    colorcet('D1')
end

% change reference
% nose channel 3
noz = bsxfun(@minus,bilat, bilat(3,:)); 
% M1 M2 = mastoids 343 344
mastoid = bsxfun(@minus,bilat, mean(bilat([343 344],:))); 
% A1 A2 = earlobes 345 346
earlob = bsxfun(@minus,bilat, mean(bilat([345 346],:))); 
reref = {bilat; noz; mastoid ;earlob};
mm = round(max(max(abs(cell2mat(reref(:))))),-1);
figure('position', [200, 1000, 2000, 1000])
for reference = 1:4
    for roi=1:9
        subplot(4,9,roi+9*(reference-1));hold on
        plotTopo(reref{reference}(:,roi),'layout/Standard-10-5-Cap385.sfp')
        caxis([-mm mm]);title(listBilat{roi}(1:end-2));
        colorcet('D1')
    end
end
saveas(gcf,['figures/bilatRef-fieldtrip'],'png')
saveas(gcf,['figures/bilatRef-fieldtrip'],'fig')
figure('position', [200, 1000, 2000, 1000])
for reference = 1:4
    for roi=1:9
        subplot(4,9,roi+9*(reference-1));hold on
        topoplot(reref{reference}(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','off' );
        caxis([-mm mm]);title(listBilat{roi}(1:end-2));
        colorcet('D1')
    end
end
saveas(gcf,['figures/bilatRef-eeglab'],'png')
saveas(gcf,['figures/bilatRef-eeglab'],'fig')


rerefNeg = {bilat.*-1; noz.*-1; mastoid.*-1 ;earlob.*-1};
figure('position', [200, 1000, 2000, 1000])
for reference = 1:4
    for roi=1:9
        subplot(4,9,roi+9*(reference-1));hold on
        topoplot(rerefNeg{reference}(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','off' );
        caxis([-mm mm]);title(listBilat{roi}(1:end-2));
        colorcet('D1')
    end
end

% absolute does not make sense: head should be 0
rerefAbs = {abs(bilat); abs(noz); abs(mastoid) ;abs(earlob)};
figure('position', [200, 1000, 2000, 1000])
for reference = 1:4
    for roi=1:9
        subplot(4,9,roi+9*(reference-1));hold on
        topoplot(rerefAbs{reference}(:,roi),'layout/Standard-10-5-Cap385.sfp','colormap',colorcet('D1'),'electrodes','off' );
        caxis([-mm mm]);title(listBilat{roi}(1:end-2));
        colorcet('D1')
    end
end
    
    