function [sourceAmplitude, sourceSSVEP, sourceERP,timeERP] = createSourceROI(numROIs,diffSource,sameSource)
% for each diffSource, generate a random amplitude activity
% (sourceAmplitude) and a timecourse over 90 points (sourceERP) 
% + timeERP = time window of the ERP (2 std)

% diffSource is the indexes of the sources to generate (e.g. in the left hemisphere)
% sameSource [optional] is the indexes of the same sources as in diffSources (e.g. in the right hemisphere)
% For example, say V1L is index 1 and V1R is index 2, diffSource will be =1 and sameSource = 2. 
% The 2 vectors must respect the order of the sources and have the same length
% numROIs is the total number of ROIs

% time line
x = 0 : pi / 45 : 2 * pi-pi/45; % 360 deg with point every 4 deg
% initialise variable
sourceAmplitude = zeros(numROIs,1);
sourceERP = zeros(numROIs, length(x) );
sourceSSVEP = zeros(numROIs, length(x) );

if ~isempty(sameSource)
    if length(diffSource)~=length(sameSource)
        fprintf('ERROR: incompatible nb of sources\n')
    end
end

for k = 1 : length(diffSource)
    % random amplitude between 1 and 10
    sourceAmplitude(diffSource(k)) = 1 + randi(90)/10;
    % create SSVEP waveform with random phase delay and random
    % coef for each source
    y = rand * cos( 2 * x - rand * pi ) + ...
        rand * cos( 4 * x - rand * pi ) + ...
        rand * cos( 8 * x - rand * pi ) + ...
        rand * cos( 12 * x - rand * pi );
    % create ERP waveform
%     yERPgauss = y.*gaussmf(x,[pi/4 pi]); % std=pi/4, mean=pi
    % use a half cosine
    cosFilt = cos(-pi:pi/45:pi-pi/45);
    cosFilt(cosFilt<0) = 0;
    yERP = y.* cosFilt;
    % multiply waveform by the source amplitude (ROIxtime)
    sourceSSVEP( diffSource(k) , : )  = y * sourceAmplitude(diffSource(k));
    sourceERP( diffSource(k) , : )  = yERP * sourceAmplitude(diffSource(k));
    if ~isempty(sameSource)
        sourceSSVEP( sameSource(k) , : )  = y * sourceAmplitude(diffSource(k));
        sourceERP( sameSource(k) , : )  = yERP * sourceAmplitude(diffSource(k));
        sourceAmplitude( sameSource(k) , : )  = sourceAmplitude(diffSource(k));
    end
end

% get timepoints corresponding to the active ERP window = 2std
% centre gaussian +/- 2 std in the gaussmf
% timeERPgauss = find(x>pi-2*pi/4 & x<pi+2*pi/4);
timeERP = find(cosFilt>0);

