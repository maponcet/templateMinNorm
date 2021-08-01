
addpath([pwd filesep 'subfunctions' filesep]);


for rep = 1:100
    
SNRlevel = 20; % 0.1 means 10 times more noise than signal
signal = repmat([zeros(1,10) randsample(-21:2:21,10,1)],6,1).*rand(6,20); % electrodes x time points
timeWin = 11:20;
baseWin = 1:10;
[noise] = add_ERPnoise_with_SNR( signal, SNRlevel, timeWin);
simul = signal+noise;

% figure;imagesc(signal);colorbar;
% figure;imagesc(noisy_data);colorbar;
% figure;imagesc(simul);colorbar;

% aa(rep) = (rms(simul)/rms(noise))^2 - 1;
% bb(rep)=(rms(simul(:,timeWin))/rms(simul(:,baseWin)))^2 -1;
% cc(rep)=(sum(sum(signal.^2+noise.^2)))/sum(sum(noise.^2)) - 1;

for elec=1:6
m1(elec,rep) = (rms(simul(elec,:))/rms(noise(elec,:)))^2 - 1;
m2(elec,rep) = (rms(simul(elec,timeWin))/rms(simul(elec,baseWin)))^2 -1;
m3(elec,rep) = (sum(sum(signal(elec,:).^2+noise(elec,:).^2)))/sum(sum(noise(elec,:).^2)) - 1;
m4(elec,rep) = (rms(simul(elec,timeWin))/rms(noise(elec,timeWin)))^2 - 1;
m5(elec,rep) = (rms(signal(elec,:))/rms(noise(elec,:)))^2;
end

end

mean(mean(m1,2))
mean(mean(m2,2))
mean(mean(m3,2))
mean(mean(m4,2))
mean(mean(m5,2))

% SNR half than expected if I feed signal time window and calculate the
% output SNR on the full window.