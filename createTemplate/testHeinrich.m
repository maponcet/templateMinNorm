
load('/Users/marleneponcet/Documents/Git/templateMinNorm/createTemplate/templates/template_V1Split_Standard_1005.mat')
% easycap active electrodes 32
elecNames = {'Fp1' 'Fp2' 'F7' 'F3' 'Fz' 'F4' 'F8' 'FT10' 'FC6' 'FC2' 'FC1'...
    'FC5' 'FT9' 'T7' 'C3' 'Cz' 'C4' 'T8' 'TP10' 'CP6' 'CP2' 'CP1' 'CP5' 'TP9' 'P7' 'P3' 'Pz'...
    'P4' 'P8' 'O2' 'Oz' 'O1'}; % ref 'FCz' 
match32elec = cell2mat(arrayfun(@(x) cellfind(templates.label,elecNames{x}),1:length(elecNames),'uni',false));
allLoc = ft_read_sens('layout/Standard-10-5-Cap385.sfp');
coord = allLoc.chanpos(match32elec,:);

fid = fopen('layout/testV1Heinrich.sfp','w');
for elec=1:length(elecNames)
    fprintf(fid,'%s\t%f\t%f\t%f\t\n',elecNames{elec},coord(elec,:));
end

maxInd = round(max(abs([sum(templates.weights(match32elec,[1 2]),2); sum(templates.weights(match32elec,[3 4]),2)])));

figure;
subplot(1,2,1)
topoplot(sum(templates.weights(match32elec,[1 2]),2),'layout/testV1Heinrich.sfp','colormap',colorcet('D1'),'electrodes','ptslabels' );
caxis([-maxInd maxInd]);
title('V1 ventral (UVF)')
subplot(1,2,2)
topoplot(sum(templates.weights(match32elec,[3 4]),2),'layout/testV1Heinrich.sfp','colormap',colorcet('D1'),'electrodes','ptslabels' );
caxis([-maxInd maxInd]);
title('V1 dorsal (LVF)')
set(gcf,'position',[200 500 700 300])
saveas(gcf,['figures' filesep 'testV1Heinrich'],'png')
