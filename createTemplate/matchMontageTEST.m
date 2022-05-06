function customTemplates = matchTemplate(channelInfo,refIndex,varargin)
% function creating ROI templates corresponding to the EEG montage with 
% which the data was recorded

% input channelInfo: variable containing electrode definition
% refIndex: reference of the EEG montage. Enter the channel number or 0 
% for average reference

% For EEGlab users = EEG.chanlocs (obtained from readlocs or using the GUI)
% https://eeglab.org/tutorials/04_Import/Channel_Locations.html 
% For Fieldtrip users, structure (eg obtained from ft_read_sens) containing
%   elec.elecpos = Nx3 matrix with carthesian (x,y,z) coordinates of each
%                  electrode
%   elec.label   = cell-array of length N with the label of each electrode
%   elec.chanpos = Nx3 matrix with coordinates of each sensor

% To do the fitting, need 3-D cartesian coordinates of the electrodes
% corresponding to the EEG montage that is used AND at least 8 channel 
% labels to match. Can be from the dataset or entered manually when prompted. 

% The ROI templates were created based on 3D electrode coordinates from the 
% fieldtrip standard_1005.elc file constructed by Robert Oostenveld.
% The file can be read with ft_read_sens. 
% The electrode positions are represented in mm in the MNI coordinate 
% system: RAS for first dimension orients towards Right, the 
% second dimension orients towards Anterior, the third dimension orients 
% towards Superior. 
% see also: https://www.fieldtriptoolbox.org/template/electrode/
% for coordinate system: https://www.fieldtriptoolbox.org/faq/coordsys/
% EEGlab uses a different coordinate system: ALS. x is towards the nose, 
% y is towards the left ear, and z towards the vertex. So before finding 
% the best alignment between the dataset and templates, the coordinates 
% need to be flipped accordingly. 
% The coordinate system is by default ALS for EEGlab users, RAS otherwise. 
% 'ALS' or 'RAS' can be specified as an optional input if required (other 
% coordinate system can potentially be added in the program) 

% template assumed to be in the same directory, folder template but
% directory can be picked manually when prompted if cannot be found


if length(channelInfo) == 1 % not eeglab
    eeglabUser = 0; coordsys = 'RAS';nbChan = length(channelInfo);
else
    eeglabUser = 1; coordsys = 'ALS';nbChan = length(channelInfo.label);
end
if varargin == 1 % use coordinate system entered by user
    coordsys = varargin{1};
end
fprintf('Assumes %s coordinate system.\n',coordsys)

% try to match electrodes between dataset and template using standard 10-20: Fp1 Fp2 Fz F7 F3 C3 T7 P3 
% P7 Pz O1 Oz O2 P4 P8 T8 C4 F4 F8 Cz
listLabels = {'Fp1' 'Fp2' 'Fz' 'F7' 'F3' 'C3' 'T7' 'P3' 'P7' 'Pz' 'O1' 'Oz' ...
    'O2' 'P4' 'P8' 'T8' 'C4' 'F4' 'F8' 'Cz' };

% get chan indexes that match between the listLabels and the data
if eeglabUser
    matchIndex = arrayfun(@(x) cellfind({channelInfo.labels},listLabels{x}),1:length(listLabels),'uni',false);
else
    matchIndex = arrayfun(@(x) cellfind(channelInfo.label,listLabels{x}),1:length(listLabels),'uni',false);
end

nbMatch = sum(~cellfun('isempty', matchIndex)); % find nb of matching labels
% emptyIndex = cellfun('isempty', matchIndex); % find indexes of labels not present in the data
% matchIndex(emptyIndex) = {0}; % set absent channels to 0
matchIndex = cell2mat(matchIndex); % convert cell to matrix
    
% if data does not contain channel labels
% ask user to enter a set of channel labels with their corresponding
% indexes from standard 10-20
while nbMatch == 0
    prompt = {'Enter the channel numbers for the following channels in the SAME order, separated by space. If a channel is not in your EEG montage write 0. '...
        'Fp1 Fp2 Fz F7 F3 C3 T7 P3 P7 Pz O1 Oz O2 P4 P8 T8 C4 F4 F8 Cz'};
    dlgtitle = 'Create ROI-templates corresponding to the dataset';
    answer = inputdlg(prompt,dlgtitle);
    chanIndex = str2num(answer{1});
    if length(chanIndex)~=length(listLabels)
        msgError = errordlg('Please enter as many numbers as the number of channels in the list','Input error','modal');
        uiwait(msgError);
    else
        if sum(~ismember(chanIndex,0:label))>0 % include 0 for non-existing channel
            msgError = errordlg('Please enter valid channel number','Input error','modal');
            uiwait(msgError);
        else
            nbMatch = sum(chanIndex>0);
        end
    end
end
% Biosemi128: 93 80 85 103 0 115 119 0 127 19 15 23 28 0 43 58 54 0 71 1

% this is an arbitrary threshold for min nb of matching labels
if nbMatch >0 && nbMatch < 8 
    error('Could not find enough matching channels')
end
    
% reference of the EEG montage
while length(refIndex)~=1 || ~ismember(refIndex,0:nbChan) % include 0 for average reference
    msgError = errordlg('Please enter a valid number for the EEG reference','Input error','modal');
    uiwait(msgError);
    prompt = {'Enter the reference of the montage (enter the channel number or 0 for average reference)'};
    dlgtitle = 'Create ROI-templates corresponding to the dataset';
    answer = inputdlg(prompt,dlgtitle);
    refIndex = str2num(answer{1});
end
    
% load templates
templateFile = ['temples' filesep 'template_Standard_1005.mat'];
while ~exist(templateFile)
    uiwait(msgbox('Could not find template_Standard_1005.mat file, please select it','Information','modal'));
    templateDir = uigetdir('','Pick directory containing the templates');
    templateFile = [templateDir filesep 'template_Standard_1005.mat'];
end
load(templateFile)

% get the indexes that match with the data (& from the list of channels)
listLabelMatch = {channelInfo(matchIndex).labels};
matchIndexROI = cell2mat(arrayfun(@(x) cellfind(elecDef.label,listLabelMatch{x}),1:length(listLabelMatch),'uni',false));

% get electrodes coordinate to match to the templates
% change orientation if necessary
if eeglabUser
    if strcmp(coordsys,'ALS')
        coordToMatch = [-[channelInfo.Y]; [channelInfo.X]; [channelInfo.Z]]';
    else % RAS
        coordToMatch = [[channelInfo.X]; [channelInfo.Y]; [channelInfo.Z]]';
    end
else % fieldtrip
    if strcmp(coordsys,'ALS')
        coordToMatch = [-channelInfo.chanpos(:,2) channelInfo.chanpos(:,1) channelInfo.chanpos(:,3)];
    else
        coordToMatch = channelInfo.chanpos;
    end
end


% align the montages using the matching electrodes
affineMtx = coordToMatch(matchIndex,:)' / elecDef.chanpos(matchIndexROI,:)' ;
% Apply transform to all the template electrodes
elecTrans = (affineMtx*elecDef.chanpos')';

% Use knnsearch to find the best match between electrodes of the current 
% montage and the templates.
% returns electrode indexes (same length as nbChan)
bestElec = knnsearch(elecTrans,coordToMatch);

% select only the best matching electrodes from the templates
matchedTemplate = avMap(bestElec,:); 

% re-reference the templates
if refIndex == 0
    % average reference: substract average across channels
    customTemplates = bsxfun(@minus,matchedTemplate, mean(matchedTemplate)); 
else
    % substract the reference channel
    customTemplates = bsxfun(@minus,matchedTemplate, matchedTemplate(refIndex,:)); 
end





%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% stuff while testing - to be deleted


[ALLEEG EEG CURRENTSET ALLCOM] = eeglab;
EEG = pop_loadset('filename','VP01_rej.set','filepath','/Users/marleneponcet/Desktop/testEEGlabData/');
[ALLEEG, EEG, CURRENTSET] = eeg_store( ALLEEG, EEG, 0 );



testField = ft_read_sens('EEG systems/standard/standard_1020.elc');
elecDef = ft_read_sens('EEG systems/standard/standard_1005.elc');

testField.chanpos = testField.chanpos .* [0.9 1.2 0.95];

% plot
figure;scatter3(elecDef.chanpos(:,1),elecDef.chanpos(:,2),elecDef.chanpos(:,3),'filled')
hold on;scatter3(coordToMatch(:,1),coordToMatch(:,2),coordToMatch(:,3),'filled')
hold on;scatter3(testField2(:,1),testField2(:,2),testField2(:,3),'filled')
hold on;scatter3(elocTrans(:,1),elocTrans(:,2),elocTrans(:,3),'filled')

figure; scatter3(coordToMatch(matchIndex,1),coordToMatch(matchIndex,2),coordToMatch(matchIndex,3),'filled')
hold on;scatter3([channelInfo(matchIndex).X],[channelInfo(matchIndex).Y],[channelInfo(matchIndex).Z],'filled')
hold on;scatter3(elecDef.chanpos(matchIndexROI,1),elecDef.chanpos(matchIndexROI,2),elecDef.chanpos(matchIndexROI,3),'filled')

figure; scatter3(coordToMatch(matchIndex,1),coordToMatch(matchIndex,2),coordToMatch(matchIndex,3),'filled')
hold on; scatter3(elecTrans(:,1),elecTrans(:,2),elecTrans(:,3),'filled')
hold on; scatter3(elecTrans(matchIndex,1),elecTrans(matchIndex,2),elecTrans(matchIndex,3),'filled')

elecDef.label(bestElec)

loc = [1:9;10:18]; loc = loc(:);
mm = round(max(max(abs(customTemplates.data))),-1);
 figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        plotTopo(customTemplates.data(:,roi),'EEG systems/standard/standard_1020.elc');
        caxis([-mm mm]);
        title(listROIs(roi));
        colorcet('D1') 
    end
    
    mm = round(max(max(abs(matchedTemplate))),-1);
 figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        plotTopo(matchedTemplate(:,roi),'EEG systems/standard/standard_1020.elc');
        caxis([-mm mm]);
        title(listROIs(roi));
        colorcet('D1') 
    end

    mm = round(max(max(abs(matchedTemplate))),-1);
      figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        topoplot(matchedTemplate(:,roi),EEG.chanlocs,'colormap',colorcet('D1'),'electrodes','on' );
        caxis([-mm mm]);title(listROIs(roi));
    end
    
    mm = round(max(max(abs(customTemplates.data))),-1);
      figure('position', [200, 1000, 2000, 500])
    for roi=1:18
        subplot(2,9,loc(roi))
        topoplot(customTemplates.data(:,roi),EEG.chanlocs,'colormap',colorcet('D1'),'electrodes','on' );
        caxis([-mm mm]);title(listROIs(roi));
    end
    



