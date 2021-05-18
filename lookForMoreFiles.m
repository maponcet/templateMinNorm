
% list of the current sbj (skeri) files
currentList = dir(fullfile('/Volumes/Amrutam/Marlene/JUSTIN/topoMapForBiosemi/mrcProj', '**','skeri*'));
currentSbj=[];curr = 1;
for ff=1:length(currentList)
    if length(currentList(ff).name) == 9
        currentSbj(curr) = str2num(currentList(ff).name(end-3:end));
        curr = curr+1;
    end
end


% look in the 4D2 folder
directoryName = '/Volumes/Karjikai/4D2 part I/';
fileName = 'elp2mri.tran';

fileList = dir(fullfile(directoryName, '**', fileName));
fprintf('found %d files',length(fileList))
% folderList = dir(fullfile(directoryName, '**','skeri*'));
errorList = [];
for ff=1:length(fileList)
    subjNb = str2num(fileList(ff).folder(end-9:end-6));
    if isempty(subjNb)
        fprintf('error %d \n',ff)
        errorList = [errorList ff];
    elseif ismember(subjNb,currentSbj)==0
        fprintf('new skeri found! %d \n',subjNb)
        newSkeri = ['skeri' fileList(ff).folder(end-9:end-6)];
        mkdir(newSkeri)
        copyfile(fileList(ff).folder(1:end-5), newSkeri)
        currentSbj(end+1) = subjNb;
    end
end


% add the ones from NewSkeri
addList = dir(fullfile('/Volumes/Amrutam/Marlene/JUSTIN/NewSkeri', 'skeri*'));
for ff=1:length(addList)
    currentSbj(curr) = str2num(addList(ff).name(end-3:end));
    curr = curr+1;
end


% look in the 4D2 folder
directoryName = '/Volumes/NeepNTatty/4D2/part2/';
fileName = 'elp2mri.tran';

fileList = dir(fullfile(directoryName, '**', fileName));
fprintf('found %d files',length(fileList))
% folderList = dir(fullfile(directoryName, '**','skeri*'));
errorList = [];
for ff=1:length(fileList)
    subjNb = str2num(fileList(ff).folder(end-9:end-6));
    if isempty(subjNb)
        fprintf('error %d \n',ff)
        errorList = [errorList ff];
    elseif ismember(subjNb,currentSbj)==0
        fprintf('new skeri found! %d \n',subjNb)
        newSkeri = ['skeri' fileList(ff).folder(end-9:end-6)];
        mkdir(newSkeri)
        copyfile(fileList(ff).folder(1:end-5), newSkeri)
        currentSbj(end+1) = subjNb;
    end
end





% % copy some skeri files for biosemi model -> inv fwd Axx tran
% % USELESS! I need almost everything!! 
% listSkeri = dir('/Volumes/Amrutam/Marlene/JUSTIN/OriginalSkeriFolders/skeri*');
% dirFolder = '/Volumes/Amrutam/Marlene/JUSTIN/skeriBiosemi/';
% for ff=1:length(listSkeri)
%     mkdir(dirFolder,listSkeri(ff).name)
%     mkdir([dirFolder listSkeri(ff).name],'_MNE_')
%     mkdir([dirFolder listSkeri(ff).name],'Exp_MATL_HCN_128_Avg')
%     mkdir([dirFolder listSkeri(ff).name],'_dev_')
%     copyfile([listSkeri(ff).folder '/' listSkeri(ff).name '/_MNE_'],[dirFolder listSkeri(ff).name '/_MNE_/'])
%     copyfile([listSkeri(ff).folder '/' listSkeri(ff).name '/_dev_'],[dirFolder listSkeri(ff).name '/_dev_/'])
% end


% same Axx_c001.mat for all sbj
dirFolder = '/Volumes/Amrutam/Marlene/JUSTIN/OriginalSkeriFolders/';
listSkeri = dir([dirFolder 'skeri*']);
for ff=1:length(listSkeri)
    % delete folder
    rmdir([dirFolder listSkeri(ff).name filesep 'Exp_MATL_HCN_128_Avg'],'s')
    % create directory 
    mkdir([dirFolder listSkeri(ff).name],'Exp_MATL_HCN_128_Avg')
    % copy Axx
    copyfile('/Volumes/Amrutam/Marlene/JUSTIN/Axx_c001.mat',[dirFolder listSkeri(ff).name filesep 'Exp_MATL_HCN_128_Avg/'])
end


% look in the 4D2 folder
directoryName = '/Volumes/NeepNTatty/4D2/';
fileName = 'skeri0062';
fileList = dir(fullfile(directoryName, '**', fileName));


