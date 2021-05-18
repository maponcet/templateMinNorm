% function [ output ] = makeManyForwards(projectInfo,optns)
% create .mat files with fwdMatrix and roiInfo from the fif files

addpath (genpath('/Volumes/Amrutam/Marlene/Git/svndl_code'))
addpath /Volumes/Amrutam/Marlene/Git/fieldtrip-aleslab-fork/external/mne

sbjList = dir('/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/OriginalSkeriFolders/skeri*');

for ss=1:length(sbjList)
    
    srcSpace = readDefaultSourceSpace(sbjList(ss).name);
    
    fwdFile = fullfile([sbjList(ss).folder filesep sbjList(ss).name filesep],'_MNE_',[sbjList(ss).name '-fwd.fif']);
    fwd = mne_read_forward_solution(fwdFile);
    
    [fwdMatrix] = makeForwardMatrixFromMne(fwd,srcSpace);
    
    anatDir = getpref('mrCurrent','AnatomyFolder');
    roiDir  = fullfile(anatDir,sbjList(ss).name,'Standard','meshes','ROIs');
    roiInfo = getRoisByType(roiDir,'func');
    
    outputFilename=fullfile('forwardAllEGI',['forwardAndRois-' sbjList(ss).name '.mat']);
    save(outputFilename,'fwdMatrix','roiInfo');
    
end
