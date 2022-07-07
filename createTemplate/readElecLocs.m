function [fifIdx, elecIdx, elecCoord] = readElecLocs(eegSystem)

fifIdx = [ 22 9 11 33 24 36 45 52 58 62 70 75 83 92 96 108 104 124 122];
% matching electrodes between systems
matchElec = {'Fp1' 'Fp2' 'Fz' 'F7' 'F3' 'C3' 'T7' 'P3' 'P7' 'Pz' 'O1' 'Oz' 'O2' 'P4' 'P8' 'T8' 'C4' 'F4' 'F8'};
if ~exist('ft_read_sens','file')
    addpath (genpath('/Users/marleneponcet/Documents/Git/fieldtrip-aleslab-fork/'))
end

switch eegSystem
     
    case 'Standard_1005' 
        A = ft_read_sens('templates/standard_1005.elc');
        elecCoord = A.elecpos;
        elecIdx = zeros(length(matchElec),1);
        for mm = 1:length(matchElec)
            elecIdx(mm) = find(strcmp(A.label,matchElec(mm)));
        end
        
    case 'EGI32' 
        elecIdx= [1 2 17 11 3 5 13 7 17 19 9 20 10 8 16 14 6 4 12];        
        A = ft_read_sens('templates/GSN-HydroCel-32.sfp');
        elecCoord = A.elecpos;
        
%         fid = fopen('EEG systems/EGI/AdultAverageNet32_v1.sfp');
%         A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
%         A = cell2mat(A);
%         elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
%         fclose(fid);  
        
    case 'EGI64' 
        elecIdx= [ 10 5 6 18 12 20 24 28 30 34 35 37 39 42 44 52 50 60 58];
        A = ft_read_sens('templates/GSN-HydroCel-64.sfp');
        elecCoord = A.elecpos;
        
%         fid = fopen('EEG systems/EGI/AdultAverageNet64_v1.sfp');
%         A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
%         A = cell2mat(A);
%         elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
%         fclose(fid);
        
    case 'EGI128' 
        elecIdx= fifIdx;
        A = ft_read_sens('templates/GSN-HydroCel-128.sfp');
        elecCoord = A.elecpos;
        %         fid = fopen('EEG systems/EGI/GSN-HydroCel-129.sfp');
%         fid = fopen('EEG systems/EGI/AdultAverageNet128_v1.sfp');
%         A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
%         A = cell2mat(A);
%         elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
%         fclose(fid);
    
    case 'EGI256'
        elecIdx= [37 18 21 47 36 59 69 87 96 101 116 126 150 153 170 202 183 224 2];
        A = ft_read_sens('templates/GSN-HydroCel-256.sfp');
        elecCoord = A.elecpos;
        
%         fid = fopen('EEG systems/EGI/AdultAverageNet256_v1.sfp');
%         A = textscan(fid,'%*s %f %f %f'); % skip 1st column of names
%         A = cell2mat(A);
%         elecCoord = A(4:end-1,:); % skip 1st 3 fiducials & Cz
%         fclose(fid);
        
                
    otherwise
        disp('Sorry montage not currently available!')
end
