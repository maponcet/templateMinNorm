function [snr,snrPrct] = computeMetricsSNR_test(beta,source,winBaseline)

if length(winBaseline) == length(beta)/2
    winSig = max(winBaseline)+1:length(beta);
    
%     % compute SNR for each ROI
%     % SNR = (winSig - winBaseline) / winBaseline = winSig/ winBaseline - 1
%     % since winSig contains both signal + noise
%     energySNR = zeros(size(beta,1),1);
%     for roi=1:size(beta,1)
%         energySNR(roi) = (rms(beta(roi,winSig)) / rms(beta(roi,winBaseline)) )^2 - 1;
%     end
%     activ_sources = source(:,winSig(1))~=0; % ttime=1 since I only need one vector (the rest of the window is/should be the same)
%     pctEnergy = sum(energySNR(activ_sources)) / sum(energySNR);
%     %     snr = mean(mean(energySNR));
%     
%     temp1 = (rms(beta(activ_sources,winSig)) / rms(beta(activ_sources,winBaseline))) ^2 - 1;
%     temp2 = (rms(beta(:,winSig)) / rms(beta(:,winBaseline))) ^2 - 1;
%     pctEnergy2 = temp1 / temp2;
%     
%     srcPct = sum(beta(activ_sources,winSig)) / (sum(beta(activ_sources,winSig)) + sum(beta(activ_sources,winBaseline)));
%     nosources = source(:,winSig(1))==0;
%     srcPctNosrc = sum(beta(nosources,winSig)) / (sum(beta(nosources,winSig)) + sum(beta(nosources,winBaseline)));

        activ_sources = source(:,winSig(1))~=0;
        sigNoise = beta(activ_sources,winSig);
        noise = beta(activ_sources,winBaseline);
        all = sigNoise + noise;
        % A/C
        snr = (rms(sigNoise(:)) / rms(noise(:)) )^2 -1;
        % A/(A+C)
        snrPrct = (rms(sigNoise(:)) / rms(all(:)) )^2 ;

% activ_sources = source(:,winSig(1))~=0;
% normBeta = abs(beta(:,11)) / max(beta(:,11));
% srcPct = sum(beta(activ_sources,winSig)) / (sum(beta(activ_sources,winSig)) + sum(beta(activ_sources,winBaseline)));
% sum(beta(activ_sources,winSig)) / sum(beta(activ_sources,winBaseline))
% 
% norm_beta(:,nT) = beta(:,nT) / max( abs(beta(:,nT)) ); % normalise estimated sources
% relEnergy(nT) = sum( abs( norm_beta(activ_sources,nT) ) ) / sum( abs(norm_beta(:,nT)) );

else
    %%%%%%%% if more than 1 signal window
    % signal is V1, MT, V1+MT
    % each of these 3 time windows are equal to 1 baseline time window
    % compute SNR separately for the 3 signal time windows
    
    sizeWin = length(winBaseline);
    nbWinSignal = length(beta) / sizeWin - 1; % total windows of baseline size (remove 1 for only signal window);
    winSig = zeros(nbWinSignal,sizeWin);
    for ww=1: nbWinSignal
        winSig(ww,:) = (sizeWin*ww+1) : (sizeWin*(ww+1));
    end
   
    % should be able somewhat to check that the snr level matches between the
    % simulated and the computed one here
%     for ww=1:nbWinSignal
%         activ_sources = source(:,winSig(ww,1))~=0; % ttime=1 since I only need one vector (the rest of the window is/should be the same)
%         tmpSrcLevel(ww,:) = sum(abs(beta(activ_sources,winSig(ww,:)))) ;
%         tmpSrc(ww,:) = sum(abs(beta(activ_sources,winSig(ww,:)))) / (sum(abs(beta(activ_sources,winSig(ww,:))))...
%             + sum(abs(beta(activ_sources,winBaseline))) );
%   tmpSrc2(ww) = sum(sum(abs(beta(activ_sources,winSig(ww,:))))) / (sum(sum(abs(beta(activ_sources,winSig(ww,:)))))...
%             + sum(sum(abs(beta(activ_sources,winBaseline)))) );
%         nosources = source(:,winSig(ww,1))==0;
%         tmpNoSrc(ww,:) = sum(abs(beta(nosources,winSig(ww,:)))) / (sum(abs(beta(nosources,winSig(ww,:)))) ...
%             + sum(abs(beta(nosources,winBaseline))));
%     end
%     srcLevel = mean(tmpSrcLevel);
%     srcPct = mean(tmpSrc);
%     srcPctNosrc = mean(tmpNoSrc);

    
    for ww=1:nbWinSignal
        activ_sources = source(:,winSig(ww,1))~=0;
        sigNoise = beta(activ_sources,winSig(ww,:));
        noise = beta(activ_sources,winBaseline);
        all = sigNoise + noise;
        % A/C
        snrWin(ww) = (rms(sigNoise(:)) / rms(noise(:)) )^2 -1;
        % A/(A+C)
        snrPrctWin(ww) = (rms(sigNoise(:)) / rms(all(:)) )^2 ;
    end
    snr = mean(snrWin);
    snrPrct = mean(snrPrctWin);
    
    % C should be around A/planned SNR
    
    
%     % signal / (baseline)
%     for ww=1:nbWinSignal
%         % compute SNR for each ROI
%         % SNR = (winSig - winBaseline) / winBaseline = winSig/ winBaseline - 1
%         % since winSig contains both signal + noise
%         for roi=1:size(beta,1)
%             energySNR(roi,ww) = (rms(beta(roi,winSig(ww,:))) / rms(beta(roi,winBaseline)) )^2 - 1;
%         end
%         activ_sources = source(:,winSig(ww,1))~=0; % ttime=1 since I only need one vector (the rest of the window is/should be the same)
%         relEnergy(ww) = sum(energySNR(activ_sources,ww)) / sum(energySNR(:,ww));
%     end
% %     snr = mean(mean(energySNR));
%     % average across the 3 signal time windows
%     pctEnergy = mean(relEnergy);
%     pctEnergy2 = pctEnergy;
end