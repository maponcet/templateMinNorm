function alignToMrcElec(montage,varargin)
% align the MRI (fif file) with the specified montage
% OVERWRITES Axx_c001.fiff with the biosemi electrode locations
% e.g. alignToMrcElec('biosemi128')

if ~isempty(varargin)
    plotElec = varargin{1}; % 1 for plotting the electrodes fitting
else
    plotElec = 0;
end
projectDir = ['/Users/marleneponcet/Documents/data/skeriDATA/skeri' montage '/'];
subjectList = dir(projectDir);

[fiffElecIdx, elecIdx, newLoc] = readElecLocs(montage);
%biosemiEloc = biosemiEloc(:,[2 1 3]); %Swap x/y to match fiff
newLoc(:,4) = 1;%Add dimension for "homogeneous coordinates"  for affine fitting.
newLocFull = newLoc;

%Take just the corresponding electrodes
newLoc = newLoc(elecIdx,:);

for iSubj = 1:length(subjectList)
    subjId = subjectList(iSubj).name;
    if subjectList(iSubj).isdir==false
        continue;
    end
    if strncmp(subjId,'.',1)
        continue;
    end
    
    fileDIR = [fullfile(projectDir,subjId,'_MNE_') filesep];
    fifFile = [fileDIR 'Axx_c001.fif'];
    if isempty(fifFile)
        error('\n Cannot find MNE file');
    end
    disp(['Processing Subject: ' subjId ])
    
    fiffEvoked = fiff_read_evoked(fifFile,1);
    datafiff = fiffEvoked;
    
    %Create a matrix of eeg locations from the fiff file
    for iElec = 1:fiffEvoked.info.nchan
        
        fiffEloc(iElec,1:3) = fiffEvoked.info.chs(iElec).eeg_loc(1:3,1)';
        fiffEloc(iElec,4) = 1; %Add dimension for "homogeneous coordinates"  for affine fitting.
    end
    
    fiffElocAll = fiffEloc;
    %Take just the corresponding electrodes
    fiffEloc = fiffEloc(fiffElecIdx,:);
    
    %Affine:   fiffChan = biosemiEloc * AffineMtx
    % affineMtx = fiffChan/biosemiEloc
    %biosemiElocTransformed =
    affineMtx = fiffEloc'/newLoc';
    
    %Now apply transform to all electrodes.
    elocTrans = (affineMtx*newLocFull')';
    
    fiffEvoked.info = rmfield(fiffEvoked.info,'chs');
    fiffEvoked.info = rmfield(fiffEvoked.info,'dig');
    fiffEvoked.info = rmfield(fiffEvoked.info,'ch_names');
    
    for iElec = 1:length(newLocFull)
        
        fiffEvoked.info.chs(iElec).eeg_loc(1:3,1) = elocTrans(iElec,1:3)';
        fiffEvoked.info.chs(iElec).eeg_loc(1:3,2) = [1 0 0]';
        fiffEvoked.info.chs(iElec).loc(1:12,1) = ...
            [elocTrans(iElec,1:3) 1 0 0 0 1 0 0 0 1]' ;
        fiffEvoked.info.dig(iElec).r(1:3) = elocTrans(iElec,1:3);
        
        % not sure what it is for but requiered in the file so just write
        % what it was but 64 lines 
        fiffEvoked.info.chs(iElec).logno = iElec;
        fiffEvoked.info.chs(iElec).kind = datafiff.info.chs(1).kind;
        fiffEvoked.info.chs(iElec).range = datafiff.info.chs(1).range;
        fiffEvoked.info.chs(iElec).cal = datafiff.info.chs(1).cal;
        fiffEvoked.info.chs(iElec).coil_type = datafiff.info.chs(1).coil_type;
        fiffEvoked.info.chs(iElec).unit = datafiff.info.chs(1).unit;
        fiffEvoked.info.chs(iElec).unit_mul = datafiff.info.chs(1).unit_mul;
        fiffEvoked.info.chs(iElec).ch_name = ['EEG ' num2str(iElec,'%03.0f')];
        fiffEvoked.info.dig(iElec).kind = datafiff.info.dig(1).kind;
        fiffEvoked.info.dig(iElec).ident = datafiff.info.dig(1).kind;
        fiffEvoked.info.ch_names{iElec} = ['EEG ' num2str(iElec,'%03.0f')];
    end
    
    fiffEvoked.info.nchan = length(newLocFull);
    fiffEvoked.evoked.epochs = zeros(length(newLocFull),length(fiffEvoked.evoked.times));
        
    if plotElec == 1
        figure;hold on;
        scatter3(elocTrans(:,1),elocTrans(:,2),elocTrans(:,3))
        scatter3(fiffEloc(:,1),fiffEloc(:,2),fiffEloc(:,3),'filled')
%         scatter3(elocTrans(1,1),elocTrans(1,2),elocTrans(1,3),'filled') %Cz biosemi
    end
    fiff_write_evoked([fileDIR 'Axx_c001.fif'],fiffEvoked);
    
end

end


