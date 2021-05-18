%%% simulation of EEG response across harmonics for different duty cycles
clearvars; close all;
ff=5; % fondamental frequency

%% Generate onset/offset response (i.e. odd and even harmonics)
nT=512; % has to be divisible by 64
wtSize = 1.2; % size of time window = 1.2 sec
t=linspace(0,wtSize,nT); 

wList= [1:1:12]'; % max Freq is 12
nF = length(wList);
wAmp = [linspace(.05,1,nF/2) linspace(1,.25,nF/2) ]'; % define the amplitude response change with frequencies (peak around 8Hz)
wAmp(2:2:end) = wAmp(2:2:end)*1;
%%%%%%%%%%%%% what is the phase for? useful for waveform, less for
%%%%%%%%%%%%% spectrum, looks better with the multiplication
wPhase = (-1*wList)+repmat([3*pi/2 4*pi/3]',nF/2,1);
wPhase = (-1*wList)+repmat([0*pi/2 0*pi/2]',nF/2,1);

wAmp = cos(wPhase).*wAmp + 1i.*sin(wPhase).*wAmp;

fourierBasis = dftmtx(nT);
invFourier  = conj(fourierBasis)/nT;
waveForm = wAmp'*fourierBasis(wList+1,:);
waveForm = real(waveForm);