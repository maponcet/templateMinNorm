% create a model for the 18 regions used in PlosOne from the original EEG
% data (not transformed for biosemi)
% This script sets up simulation paramaters and calls mrSimScript.
% Script used for mrcProj data folder

addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/commonFunctions
addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork
addpath /Volumes/Amrutam/Marlene/Git/ssvepTesting/biosemiUpdated
addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork/external/mne
addpath (genpath('/Volumes/Amrutam/Marlene/Git/svndl_code'))
addpath(genpath('/Volumes/Amrutam/Marlene/Git/vistasoft/'))
ft_defaults

% use inverse and fwd from EGI system (not biosemi)
projectDir = '/Volumes/Amrutam/Marlene/JUSTIN/OriginalSkeriFolders/';

clear params;
params = skeriDefaultSimParameters;

% roiList = getRoisByType(roiDir,'func')
params.activeRoiList = {'V1-L', 'V1-R', 'V2V-L', 'V2V-R', 'V2D-L', 'V2D-R', ...
    'V3V-L','V3V-R', 'V3D-L', 'V3D-R', 'V4-L', 'V4-R', 'V3A-L', 'V3A-R',...
    'LOC-L', 'LOC-R', 'MT-L', 'MT-R'};

%Unilateral activations (all ROIS are separate)
% params.roiHarm = num2cell(eye(length(params.activeRoiList))*1e-4);
matROI = eye(length(params.activeRoiList))*1e-4;
for mm=1:length(matROI)
    params.roiHarm{mm} = matROI(mm,1:find(matROI(mm,:)==1e-4));
end
params.stepTimeByRoi = true;
      
params.condNumber = 500;
mrSimScript(projectDir,params);

% then need to average all the maps
subjectList = dir([projectDir filesep 'skeri*']);

nsub = 1;excludeSbj=[];
for iSubj = 1:length(subjectList)
    clear Axx AxxOld;
    % load Axx
    subjId = subjectList(iSubj).name;
    fileDIR = [fullfile(projectDir,subjId,'Exp_MATL_HCN_128_Avg') filesep];
    % reformat Axx so it runs with other programs
    AxxOld = load([fileDIR 'Axx_c' num2str(500) '.mat']);
    Axx = reformatAxx(AxxOld);
    % remove 0 in the matrix
    nbElec = size(Axx.cos,1);
    tmpCos = nonzeros(getfield(Axx,'sin')); % data is in sin not cos
    AxxCos = reshape(tmpCos,nbElec,length(tmpCos)/(nbElec));
    % include in average only if all areas are present for that sbj
    if size(AxxCos,2) == length(params.activeRoiList)
        activity(:,:,nsub) = AxxCos;
        nsub = nsub+1;
    else
        excludeSbj = [excludeSbj iSubj];
    end
end
avROI.activity = mean(activity,3);
avROI.locNames = params.activeRoiList;avROI.plotDim = [3 6];
avROI.nbSbj = length(subjectList);

figure('position', [200, 0, 1500, 800])
for src = 1:size(avROI.activity,2)
    subplot(avROI.plotDim(1),avROI.plotDim(2),src)
    plotContourOnScalp(squeeze(avROI.activity(:,src)),'skeri0044','datafiles/eegdata/')
    view(20,35)
    camproj('perspective')
    axis off
    title(avROI.locNames{src})
end
saveas(gcf,'averageMap50fromSkeri','png')

% save model
save('avMapFromSkeri','avROI')
    





%%
%%%%% bootstrap
subjectList = dir([projectDir filesep 'skeri*']);
% skeri 98 missing MT and LOC, 60-67-68-72 missing LOC
keepSkeri = setdiff(1:length(subjectList),excludeSbj);

load('avMapFromSkeri.mat')
nbBoot = 1000;
model = zeros(size(avROI.activity,1),size(avROI.activity,2),length(keepSkeri));
rmseAct = zeros(nbBoot,length(keepSkeri),size(avROI.activity,2));
meanActError = zeros(nbBoot,length(keepSkeri),size(avROI.activity,2));
for maxSbj = 1:length(keepSkeri)
    for bb=1:nbBoot
        if mod(bb,100)==0
            fprintf('maxSbj%d boostrap%d \n' ,maxSbj,bb)
        end
        clear model  
        % pick randomly maxSbj sbj with replacement
        pickSS = keepSkeri(randi(length(keepSkeri),1,maxSbj));
        for iSubj = 1:length(pickSS)
            clear Axx AxxOld;
            subjId = subjectList(pickSS(iSubj)).name;
            fileDIR = [fullfile(projectDir,subjId,'Exp_MATL_HCN_128_Avg') filesep];
            % reformat Axx so it runs with other programs
            AxxOld = load([fileDIR 'Axx_c' num2str(500) '.mat']);
            Axx = reformatAxx(AxxOld);
            % remove 0 in the matrix
            nbElec = size(Axx.cos,1);
            tmpCos = nonzeros(getfield(Axx,'sin')); % data is in sin not cos
            Acos = reshape(tmpCos,nbElec,length(tmpCos)/(nbElec));
            % for averaging
            model(:,:,iSubj) = Acos;
        end
        newCos = mean(model,3);

        rmseAct(bb,maxSbj,:) = rms(newCos - avROI.activity);
        meanActError(bb,maxSbj,:) = mean(newCos - avROI.activity);
    end
end
save rmseROI.mat rmseAct avROI meanActError

figure;hold on;
set(gcf,'position',[100,100,1500,700])
for src=1:size(rmseAct,3)
    bb = prctile(rmseAct(:,:,src),97.5);
    subplot(avROI.plotDim(1),avROI.plotDim(2),src)
    plot(mean(rmseAct(:,:,src),1),'LineWidth',2)
    patch([1:size(rmseAct,2) size(rmseAct,2):-1:1],[ prctile(rmseAct(:,:,src),2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
    title(avROI.locNames{src})
    xlabel('nb sbj included')
    ylabel('rmse')
    ylim([0 0.3])
    xlim([0 50])
end
saveas(gcf,['figures' filesep 'rmseAvRoi'],'png')

%%%% plot change
figure;hold on;
set(gcf,'position',[100,100,1500,900])
for src=1:size(rmseAct,3)
    change = rmseAct(:,1:49,src) - rmseAct(:,2:50,src);
    bb = prctile(change,97.5);
    subplot(avROI.plotDim(1),avROI.plotDim(2),src)
    plot(mean(change,1),'LineWidth',2)
    title(avROI.locNames{src})
    xlabel('nb sbj included')
    ylabel('change in RMSE')
%     line([0 50],[0.001 0.001],'Color','r')
%     ylim([0 0.045])
    patch([1:size(change,2) size(change,2):-1:1],[ prctile(change,2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
    ylim([0 0.25])
    xlim([0 50])
end

%%% plot variability 
figure;hold on;
set(gcf,'position',[100,100,1500,700])
for src=1:size(rmseAct,3)
    bb = prctile(rmseAct(:,:,src),97.5) - prctile(rmseAct(:,:,src),2.5);
    subplot(avROI.plotDim(1),avROI.plotDim(2),src)
    plot(mean(rmseAct(:,:,src),1),'LineWidth',2)
    title(avROI.locNames{src})
    xlabel('nb sbj included')
    ylabel('97.5-2.5 rmse prctile')
    ylim([0 0.15])
    xlim([0 50])
end




% figure;hold on;
% for src=1:size(rmseAct,3)
%     bb = prctile(meanActError(:,:,src),97.5);
%     subplot(avROI.plotDim(1),avROI.plotDim(2),src)
%     plot(mean(meanActError(:,:,src),1),'LineWidth',2)
%     patch([1:size(meanActError,2) size(meanActError,2):-1:1],[ prctile(meanActError(:,:,src),2.5) bb(end:-1:1)],'b','FaceAlpha',.1)
%     title(avROI.locNames{src})
%     xlabel('nb sbj included')
%     ylabel('mean difference')
% end
% saveas(gcf,[modelName '-meanDiff'],'png')
