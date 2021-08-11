function makeFwd(projectName)
% create forward models and .mat files with fwdMatrix and roiInfo 
% from the Axx_c001.fif files (creates skeri0001-fwd.fif, -sph-fwd)
% e.g. makeFwd('Biosemi128')

% need addpath fieldtrip-aleslab-fork/external/mne  + genpath svndl_code??? 

% create the forward models
if ~exist(['/Users/marleneponcet/Documents/data/skeriDATA/skeri' projectName],'dir')
    error('Cannot run the function: skeri folder does not exist')
end
prepareInversesForMrc(['/Users/marleneponcet/Documents/data/skeriDATA/skeri' projectName])
 
% create .mat files
mkdir(['/Users/marleneponcet/Documents/data/skeriDATA/forward' projectName])
sbjList = dir(['/Users/marleneponcet/Documents/data/skeriDATA/skeri' projectName '/skeri*']);

for ss=1:length(sbjList)
    
    srcSpace = readDefaultSourceSpace(sbjList(ss).name);
    
    fwdFile = fullfile([sbjList(ss).folder filesep sbjList(ss).name filesep],'_MNE_',[sbjList(ss).name '-fwd.fif']);
    fwd = mne_read_forward_solution(fwdFile);
    
    [fwdMatrix] = makeForwardMatrixFromMne(fwd,srcSpace);
    
    anatDir = getpref('mrCurrent','AnatomyFolder');
    roiDir  = fullfile(anatDir,sbjList(ss).name,'Standard','meshes','ROIs');
    roiInfo = getRoisByType(roiDir,'func');
    
    outputFilename=fullfile(['/Users/marleneponcet/Documents/data/skeriDATA/forward' projectName],['forwardAndRois' projectName '-' sbjList(ss).name '.mat']);
    save(outputFilename,'fwdMatrix','roiInfo');
end

