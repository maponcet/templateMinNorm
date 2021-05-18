function [Y,source,data] = simpleSimulation(roiNames,activity,ac_sources,SNRlevel,phase)
% simulate sources of different intensity and EEG response over time

%% first the sources
Sources_distribution = zeros( 1 , length(roiNames) );
% Loop on all the activated clusters
for k = 1 : length(ac_sources)
    % Random choice of the source amplitude between 1.1 to 10... WHY??
    % in HBM paper: "important magnitude differences between the different
    % clusters"
    Sources_distribution(ac_sources(k)) = 1 + randi(90)/10;
end
source = Sources_distribution;
% source = avMap.activity .* Sources_distribution;
% sumSource = avMap.activity * Sources_distribution'; % sum over the ROIs
% plotContourOnScalp(sumSource,'skeri0044','datafiles/eegdata/')

%% then over time
x = 0 : pi / 45 : 2 * pi;
% Simul signal with phase delay
y = cos( 4 * x - phase * pi / 180 ) + ...                   
    cos( 8 * x - phase * pi / 180 ) + ...
    cos( 12 * x - phase * pi / 180 ) + ...
    cos( 16 * x - phase * pi / 180 );
data = zeros( size(Sources_distribution,2) , length(x) );
for k = 1 : length(ac_sources)
    data( ac_sources(k) , : ) = y * Sources_distribution( ac_sources(k) );
end
% Creation of the EEG signal
y_stim = activity * data;  
% add noise
noisy_data = [];
if SNRlevel > 0
    [noisy_data,noise_var] = add_noise_with_SNR( y_stim , SNRlevel ); 
    Y = y_stim + noisy_data;                           
end
