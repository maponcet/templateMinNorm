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
        percErrSet(bb,iboot) = (mean(bootDraw) - mean(setOfData)).^2 / mean(setOfData).^2; 
    end
end
% change in percent error
figure(3); clf 
loglog(mean(percErrSet,2))
hold on;
loglog(mean(sqrt(percErrSet),2))
loglog(1./sqrt(1:100))
legend('square error','RMS error','1/sqrt(N)')
title('same set: percent error')
% variability of error measure
figure(4); clf
loglog(var(percErrSet,[],2))
hold on;
loglog(std(percErrSet,[],2))
loglog(1./sqrt(1:100))
title('same set: variability')
legend('var','std','1/sqrt(N)')

figure;loglog(mean(sqrt(percErrSet),2)-mean(sqrt(percErr),2))