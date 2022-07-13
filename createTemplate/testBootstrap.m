%%%%% test bootstrap ....
%%%%
trueMean = 10;
for iboot=1:1000
    for bb=1:100
        bootDraw = 10*rand(bb,1) + trueMean;
        percErr(bb,iboot) = (mean(bootDraw) - trueMean).^2 / trueMean.^2; 
    end
end
% change in percent error
figure(1); clf 
loglog(mean(percErr,2))
hold on;
loglog(mean(sqrt(percErr),2))
legend('square error','RMS error')
title('random: percent error')
% variability of error measure
figure(2); clf
loglog(var(percErr,[],2))
hold on;
loglog(std(percErr,[],2))
loglog(1./sqrt(1:100))
legend('var','std','1/sqrt(N)')
title('random: variability')

%%%%
setOfData = rand(1,100,1)*10+5;
for iboot=1:1000
    for bb=1:100
        bootDraw = datasample(setOfData,bb);
        randomDraw = rand(1,bb)*10+5;
        percErrSet(bb,iboot) = (mean(bootDraw) - mean(setOfData)).^2 / mean(setOfData).^2; 
        percErrRand(bb,iboot) = (mean(randomDraw) - mean(setOfData)).^2 / mean(setOfData).^2; 
    end
end
% change in percent error
figure(3); clf 
loglog(mean(percErrSet,2))
hold on;
loglog(mean(sqrt(percErrSet),2))
loglog(mean(percErrRand,2))
loglog(mean(sqrt(percErrRand),2))
loglog(1./sqrt(1:100))
legend('square error','RMS error','rand square error','rand RMS error','1/sqrt(N)')
title('same set: percent error')
% variability of error measure
figure(4); clf
loglog(var(percErrSet,[],2))
hold on;
loglog(std(percErrSet,[],2))
loglog(var(percErrRand,[],2))
loglog(std(percErrRand,[],2))
loglog(1./sqrt(1:100))
title('same set: variability')
legend('var','std','rand var','rand std','1/sqrt(N)')


%%%% Different mean
% total error = bias + variance. bias is constant while variance decreases
% by increasing number of samples. If big bias (ie between 2 means) then 
% error dominated by bias, if same mean then dominated by variance. So
% adding more sbj does not make much difference after 50. There might be a
% bias due to age, ethnicity etc. but as there would be in any expt.
setOfData = rand(1,100,1)*10+5;
for iboot=1:1000
    for bb=1:100
        bootDraw = datasample(setOfData,bb);
        randomDraw = rand(1,bb)*20+5;
        percErrSet(bb,iboot) = (mean(bootDraw) - mean(setOfData)).^2 / mean(setOfData).^2; 
        percErrRand(bb,iboot) = (mean(randomDraw) - mean(setOfData)).^2 / mean(setOfData).^2; 
    end
end
% change in percent error
figure; clf 
loglog(mean(percErrSet,2))
hold on;
loglog(mean(sqrt(percErrSet),2))
loglog(mean(percErrRand,2))
loglog(mean(sqrt(percErrRand),2))
loglog(1./sqrt(1:100))
legend('square error','RMS error','rand square error','rand RMS error','1/sqrt(N)')
title('diff set: percent error')




%%%% add one dim for electrodes (3 elec)
setOfDataElec = rand(100,3)*10+5;
for iboot=1:1000
    for bb=1:100
        bootDraw = datasample(setOfDataElec,bb); % resample sbj keeping elec consistent with sbj
        randomDraw = rand(bb,3)*10+5;
        percErrSet(bb,iboot) = (mean(bootDraw) - mean(setOfDataElec)).^2 / mean(setOfDataElec).^2; 
        percErrRand(bb,iboot) = (mean(randomDraw) - mean(setOfDataElec)).^2 / mean(setOfDataElec).^2; 
    end
end
% change in percent error
figure; clf 
loglog(mean(percErrSet,2))
hold on;
loglog(mean(sqrt(percErrSet),2))
loglog(mean(percErrRand,2))
loglog(mean(sqrt(percErrRand),2))
loglog(1./sqrt(1:100))
legend('square error','RMS error','rand square error','rand RMS error','1/sqrt(N)')
title('same set: percent error')
