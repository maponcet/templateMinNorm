% compare 64 channels with locations provided by biosemi vs. from Standard
% 10-05 system
load('compareElec')
rejectElec = [6 15 41 52]; % non-standard biosemi 
plotElec = setdiff(1:64, rejectElec);
figure; hold on;
scatter3(testElecLoc.biosemi(plotElec,1), testElecLoc.biosemi(plotElec,2), testElecLoc.biosemi(plotElec,3),'filled')
scatter3(-testElecLoc.eeglab(plotElec,2), testElecLoc.eeglab(plotElec,1), testElecLoc.eeglab(plotElec,3),'filled')
legend('biosemi','standard')
saveas(gcf,'figures/compareLoc','fig')

% I have also made a computation of all the electrode locations on a 
% reallistical head surface, based on the distances along the (triangulated) 
% surface of the head. The head surface used was constructed from the 
% canonical MRI that is included in the SPM2 package, and locations are 
% expressed in MNI coordinates.
% https://robertoostenveld.nl/electrode/

% can see that for biosemi, coordinates are symetrical along the midline
% (e.g. Fp1 vs Fp2) which is not the case for the more "realistic" 
% standard 10-05 

testEEGlab = readlocs()

addpath(genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork/'));
elecEGI = ft_read_sens('EEG systems/EGI/AdultAverageNet128_v1.sfp');
elecStand = ft_read_sens('EEG systems/standard/standard_1005.elc');
elecEGIave = ft_read_sens('EEG systems/EGI/GSN-HydroCel-128.sfp');

matchIndex = arrayfun(@(x) cellfind(elecEGI.label,elecStand.label{x}),1:length(elecStand.label),'uni',false);

figure; hold on;
scatter3(elecEGI.elecpos(:,1)*10, elecEGI.elecpos(:,2)*10, elecEGI.elecpos(:,3)*10,'filled')
scatter3(elecEGIave.elecpos(:,1)*10, elecEGIave.elecpos(:,2)*10, elecEGIave.elecpos(:,3)*10,'filled')
scatter3(elecStand.elecpos(:,1), elecStand.elecpos(:,2), elecStand.elecpos(:,3),'filled')
legend('EGIaverage(used for templates)','EGI GSN','stand346')
saveas(gcf,'figures/compareEGI','fig')
