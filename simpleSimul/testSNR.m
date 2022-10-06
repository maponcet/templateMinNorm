
addpath([pwd filesep 'subfunctions' filesep]);


for rep = 1:1000
    
SNRlevel = 10; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
signal = repmat([zeros(1,10) randsample(-21:2:21,10,1)],6,1).*rand(6,20); % electrodes x time points
timeWin = 11:20;
baseWin = 1:10;
[noise] = add_ERPnoise_with_SNR( signal, SNRlevel, timeWin);
simul = signal+noise;

% figure;imagesc(signal);colorbar;
% figure;imagesc(noise);colorbar;
% figure;imagesc(simul);colorbar;

erpWin = simul(:,timeWin); preWin = simul(:,baseWin);
noiseERP = noise(:,timeWin);

aa(rep) = (rms(simul(:))/rms(noise(:)))^2 - 1;
bb(rep)= (rms(erpWin(:))/rms(preWin(:)))^2 -1;
cc(rep)= (sum(sum(signal(:).^2+noise(:).^2)))/sum(sum(noise(:).^2)) - 1;
dd(rep) = (rms(signal(:))/rms(noise(:)))^2 ;

end

mean(aa)
mean(bb)
mean(cc)
mean(dd)

% in add_ERPnoise_with_SNR, when "+ rand * variance_m" is there,  bb = SNR
% but aa & cc are half. without the randomisation part, aa = SNR and bb is
% twice SNR. 

% half of the entire window does not have signal so when using simul, the
% SNR is half. Should only use window where there is actually a signal

% noise_var(k) = 1/SNR_lin * (0.5 * variance_m + rand * variance_m)
% adds equivalent to 1*variance of the signal scales by the SNR 
% since 0.5+rand = 1 on average. 