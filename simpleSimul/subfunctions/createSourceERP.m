function sourceERP = createSourceERP(numROIs,sourceROI,sameSource)
% generate for each sourceROI a timecourse (sourceERP) over time (45ms) 
% sameSource [optional] will duplicate the generated ERP to these sources
% (indexes). Eg: sourceERP = createSourceERP(18,[1 17],[2 18]) for 18
% sources with source 1 & 17 (V1-L & MT-L) the same as 2 & 18 (V1-R & MT-R)

if nargin == 2
    sameSource = [];
end

% time line
x = 0 : pi / 45 : pi-pi/45; % 360 deg 
% initialise variable
sourceERP = zeros(numROIs, length(x) );

for k = 1 : length(sourceROI)
    % random amplitude between 1 and 10
    sourceAmplitude = 1 + randi(90)/10;
    % create SSVEP waveform with random phase delay and random
    % coef for each source
    y = rand * cos( 2 * x - rand * pi ) + ...
        rand * cos( 4 * x - rand * pi ) + ...
        rand * cos( 8 * x - rand * pi ) + ...
        rand * cos( 12 * x - rand * pi );
    % create ERP waveform
%     yERPgauss = y.*gaussmf(x,[pi/4 pi]); % std=pi/4, mean=pi
%     % use a half cosine
%     cosFilt = cos(-pi:pi/45:pi-pi/45);
%     cosFilt(cosFilt<0) = 0;
    cosFilt = cos(-pi:pi/45:pi-pi/45);
    % keep only positive nb
    filt = cosFilt(find(cosFilt>0));
    yERP = y.* filt;
    % multiply waveform by the source amplitude (ROIxtime)
    sourceERP( sourceROI(k) , : )  = yERP * sourceAmplitude;
    if ~isempty(sameSource)
        sourceERP( sameSource(k) , : )  = yERP * sourceAmplitude;
    end
end


