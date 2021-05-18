function [signal, signalERP, sourceOverTime, sourceOverTimeERP] = simulSource(areaMap,noiseLevel, diffSource,sameSource)
% Simulate a SSVEP and ERP signal over time (90 timepoints)
% return the source amplitude at the ROI for both signals
% needs the map of the roi for each sbj (matrix of electrodes * areas * sbj)
% the index of roi sources. To respect retinotopy use diffSource (create a
% new source) and sameSource(which is a copy of another source - the one at
% the same place in the vector). For example, say V1L is index 1 and V1R is 
% index 2, diffSource will be =1 and sameSource = 2. The 2 vectors must be
% of the same length. If sameSource is empty then just creates different
% signals for each of the source

x = 0 : pi / 45 : 2 * pi-pi/45;
sourceOverTime = zeros(size(areaMap,2), length(x) );
sourceOverTimeERP = zeros(size(areaMap,2), length(x) );

if ~isempty(sameSource)
    if length(diffSource)~=length(sameSource)
        fprintf('ERROR: incompatible nb of sources\n')
    end
end

for k = 1 : length(diffSource)
    % random amplitude between 1 and 10
    sourceAmplitude = 1 + randi(90)/10;
    % create SSVEP waveform with random phase delay and random
    % coef for each source
    y = rand * cos( 4 * x - rand * pi ) + ...
        rand * cos( 8 * x - rand * pi ) + ...
        rand * cos( 12 * x - rand * pi ) + ...
        rand * cos( 16 * x - rand * pi );
    % create ERP waveform (use cos envelop instead?)
    yERP = y.*gaussmf(x,[pi/5 pi]);
    % multiply waveform by the source amplitude (ROIxtime)
    sourceOverTime( diffSource(k) , : )  = y * sourceAmplitude;
    sourceOverTimeERP( diffSource(k) , : )  = yERP * sourceAmplitude;
    if ~isempty(sameSource)
        sourceOverTime( sameSource(k),:) = y * sourceAmplitude;
        sourceOverTimeERP( sameSource(k) , : )  = yERP * sourceAmplitude;
    end
end


% Simul EEG signal from sourceValOverTime for each sbj (different
% ROI map)
for iSub=1:size(areaMap,3)
    y_stim = areaMap(:,:,iSub) * sourceOverTime;
    y_stimERP = areaMap(:,:,iSub) * sourceOverTimeERP;
    % add noise
    [noisy_data,~] = add_noise_with_SNR( y_stim , noiseLevel );
    signal(128*(iSub-1)+1:128*iSub,:) = y_stim + noisy_data;
    [noisy_data,~] = add_noise_with_SNR( y_stimERP , noiseLevel );
    signalERP(128*(iSub-1)+1:128*iSub,:) = y_stimERP + noisy_data;
end
