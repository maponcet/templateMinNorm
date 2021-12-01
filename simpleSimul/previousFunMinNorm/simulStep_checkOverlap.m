clearvars;close all;
% compare variability from sbj and from noise amount
% list of sbj changes across bootstrap but same across noise levels
addpath('/Users/marleneponcet/Documents/Git/svndl_code/alesToolbox');

addpath([pwd filesep 'subfunctions' filesep]);
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
% dataPath = '/Volumes/Amrutam/Marlene/JUSTIN/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
% dirList = dirList([1:4 6:17 19:21 23:47 49 50]);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)

numROIs = length(listROIs);

% some parameters
SNRlevel = 10000; % 0.1 means 10 times more noise than signal, 10 means 10 times more signal than noise
nLambdaRidge = 10; % for calculating minimum_norm, reg constant, hyper param in min norm
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L'};
sourceR = {'V1-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

nbSbjToInclude = 25;
% nbSbjToInclude =20;

totBoot = 20;

% initialise variables
aucAve = zeros(totBoot,length(SNRlevel));
energyAve = zeros(totBoot,length(SNRlevel));
mseAve = zeros(totBoot,length(SNRlevel));
aucAveO = zeros(totBoot,length(SNRlevel));
energyAveO = zeros(totBoot,length(SNRlevel));
mseAveO = zeros(totBoot,length(SNRlevel));
aucAveNO = zeros(totBoot,length(SNRlevel));
energyAveNO = zeros(totBoot,length(SNRlevel));
mseAveNO = zeros(totBoot,length(SNRlevel));



numSubs = nbSbjToInclude;


for repBoot=1:totBoot
    fprintf('N%d bootstrap %d\n',numSubs,repBoot)
% list of random sbj with replacement
% listSub = randi(length(dirList),numSubs,1);
    listSub = randperm(length(dirList),numSubs);
    otherSbj = setdiff(1:length(dirList),listSub);

%% LOAD FWD
fullFwd=cell(1,numSubs);roiFwd=cell(numSubs,numROIs);idxROIfwd=cell(numSubs,numROIs);
for iSub=1:numSubs
    clear fwdMatrix roiInfo
    % fwd file
    load([dataPath dirList(listSub(iSub)).name])
    fullFwd{iSub} = fwdMatrix;
    %             indexROI = cell2mat(arrayfun(@(x) cellfind({roiInfo.name},listROIs{x}),1:length(listROIs),'uni',false));
    % go through each ROI and save the corresponding fwdMesh values
    % corresponding to the indexes of that ROI
    for rr=1:numROIs
        indexROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiFwd{iSub,rr} = fwdMatrix(:,roiInfo(indexROI).meshIndices);
        % to get roiFwd for one sbj= [roiFwd{iSub,:}]
        % save the index for each ROI
        idxROIfwd{iSub,rr} = roiInfo(indexROI).meshIndices;
    end
end

%% Simulate sources
% 45ms baseline + 45ms signal 15ms V1 15ms MT 15ms both
x = 0 : pi / 45 : 2 * pi-pi/45; % 360 deg with point every 4 deg
srcAmp = zeros(numROIs,1);
srcERP = zeros(numROIs, length(x) );
% V1 left & right
srcAmp(ac_sources) = 1;
srcERP([ac_sources(1) ac_sources(2)],[46:90]) = 1;
% srcERP([ac_sources(2) ac_sources(4)],[61:90]) = 1;
winERP = 46:90;
% ERP baseline timewindow
timeBase = 1:45;


for level=1:length(SNRlevel)
% fprintf('noise %d\n',level)
%%
noiseLevel = SNRlevel(level);

%%% Simulate scalp activity (Y)
% use the generated sources to simulate scalp activity for each sbj
% (using individual fwd model)
Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
Y_noise = Y;
Y_avg = zeros(numSubs,size(fullFwd{1},1),length(srcERP));

for iSub=1:numSubs
    % initialise matrix of source activity
    sourceData = zeros(size(fullFwd{iSub},2) , length(srcERP));
    sourceNoise = sourceData; % no signal, used to compute SNR
    %             if length([idxROIfwd{iSub,ac_sources}]) ~= length(unique([idxROIfwd{iSub,ac_sources}]))
    %                 fprintf('Overlapping source indexes S%d \n',listSub(iSub));
    %             end
    for ss=1:length(ac_sources)
        % note that if there is overlapping index (same idx for 2
        % ROIs), the value in sourceData will be of the latest
        % source
        sourceData(idxROIfwd{iSub,ac_sources(ss)},:) = repmat(srcERP(ac_sources(ss),:),length(idxROIfwd{iSub,ac_sources(ss)}),1);
    end
    % multiply fwd (128*20484) with the activated idx over time
    % (sourceData of 20484*90) and obtain Y elec x time
    y_stim = fullFwd{iSub} * sourceData;
    % add noise
    [noisy_data] = add_ERPnoise_with_SNR( y_stim , noiseLevel,winERP );
    % to keep the same SNR for the 2 Y, need to compute noise for
    % the 2 Y separately as it is based on the variance of the signal
    Y(iSub,:,:) = y_stim + noisy_data;
    %             Y_SSVEPprev(iSub,:,:) = y_stimSSVEP + noisy_data;
    Y_noise(iSub,:,:) = noisy_data;
    %             Y_noiseSSVEP(iSub,:,:) = noisy_dataSSVEP;
    %             Ypure(iSub,:,:) = y_stim;
    %             Y2(iSub,:,:) = y_stim + noisy_dataSSVEP;
end

%%% Use average reference for centering Y
%%% that is: substract the average electrode activity at each time point
% this is done by bsxfun which applies element-wise substraction (the 90
% averages across electrodes) - Useless
for iSub=1:numSubs
    Y_avg(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
end



%% compute minimum norm using 3 different templates
% stack all the participants vertically
stackY = reshape(permute(Y_avg,[2,1,3]),[size(Y_avg,1)*size(Y_avg,2),size(Y_avg,3)]);

% compute the template for the 25 sbj overlap 
for iSub=1:numSubs
    for rr=1:length(listROIs)
        roiMapOverlap(:,rr,iSub) = sum(roiFwd{iSub,rr},2); 
    end
end
avMapOverlap = mean(roiMapOverlap,3);

% 25 no overlap
for iSub=1:numSubs
    clear fwdMatrix roiInfo
    load([dataPath dirList(otherSbj(iSub)).name])
    for rr=1:numROIs
        clear idxROI
        idxROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiMapNoOverlap(:,rr,iSub) = sum(fwdMatrix(:,roiInfo(idxROI).meshIndices),2);
    end
end
avMapNoOverlap = mean(roiMapNoOverlap,3);

% % Compare the ROI templates = they are identical but different scales!
% figure;
% for rr=1:length(listROIs)
%     subplot(3,6,rr);plotOnEgi(avMapOverlap(:,rr));colorbar
% end
% figure;
% for rr=1:length(listROIs)
%     subplot(3,6,rr);plotOnEgi(avMapNoOverlap(:,rr));colorbar
% end
% figure;subplot(3,1,1);imagesc(avMapOverlap);colorbar;title('roi overlap')
% subplot(3,1,2);imagesc(avMapNoOverlap);colorbar;title('roi no overlap')
% subplot(3,1,3);imagesc(avMapOverlap-avMapNoOverlap);colorbar;title('difference')

% compute min norm
stackAvMap = repmat(avMap,numSubs,1);
[average(repBoot).beta, average(repBoot).betaMinNorm, average(repBoot).lambda, ...
    average(repBoot).gcvError, average(repBoot).lambdaGridMinNorm] = minimum_norm(stackAvMap, stackY, nLambdaRidge);
stackOverlap = repmat(avMapOverlap,numSubs,1);
[overlap(repBoot).beta, overlap(repBoot).betaMinNorm, overlap(repBoot).lambda, ...
    overlap(repBoot).gcvError, overlap(repBoot).lambdaGridMinNorm] = minimum_norm(stackOverlap, stackY, nLambdaRidge);
stackNoOverlap = repmat(avMapNoOverlap,numSubs,1);
[noOverlap(repBoot).beta, noOverlap(repBoot).betaMinNorm, noOverlap(repBoot).lambda, ...
    noOverlap(repBoot).gcvError, noOverlap(repBoot).lambdaGridMinNorm] = minimum_norm(stackNoOverlap, stackY, nLambdaRidge);


% stackAvMapCut = repmat(avMap(:,[1:4 13:18]),numSubs,1);
% [average(repBoot).betaCut, average(repBoot).betaMinNormCut, average(repBoot).lambdaCut, ...
%     average(repBoot).gcvErrorCut, average(repBoot).lambdaGridMinNormCut] = minimum_norm(stackAvMapCut, stackY, nLambdaRidge);
% stackOverlapCut = repmat(avMapOverlap(:,[1:4 13:18]),numSubs,1);
% [overlap(repBoot).betaCut, overlap(repBoot).betaMinNormCut, overlap(repBoot).lambdaCut, ...
%     overlap(repBoot).gcvErrorCut, overlap(repBoot).lambdaGridMinNormCut] = minimum_norm(stackOverlapCut, stackY, nLambdaRidge);
% stackNoOverlapCut = repmat(avMapNoOverlap(:,[1:4 13:18]),numSubs,1);
% [noOverlap(repBoot).betaCut, noOverlap(repBoot).betaMinNormCut, noOverlap(repBoot).lambdaCut, ...
%     noOverlap(repBoot).gcvErrorCut, noOverlap(repBoot).lambdaGridMinNormCut] = minimum_norm(stackNoOverlapCut, stackY, nLambdaRidge);

sharedVal(repBoot).avMapNoOverlap = avMapNoOverlap;
sharedVal(repBoot).avMapOverlap = avMapOverlap;
sharedVal(repBoot).stackY = stackY;
sharedVal(repBoot).roiMapOverlap = roiMapOverlap;
sharedVal(repBoot).roiMapNoOverlap = roiMapNoOverlap;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute auc, mse, relative energy using average signal in rois
% do for all the min norm outputs
[aucAve(repBoot), energyAve(repBoot), mseAve(repBoot)] = computeMetrics(average(repBoot).beta(:,winERP),srcERP(:,winERP));
[aucAveO(repBoot), energyAveO(repBoot), mseAveO(repBoot)] = computeMetrics(overlap(repBoot).beta(:,winERP),srcERP(:,winERP));
[aucAveNO(repBoot), energyAveNO(repBoot), mseAveNO(repBoot)] = computeMetrics(noOverlap(repBoot).beta(:,winERP),srcERP(:,winERP));

% [aucAve_cut(repBoot), energyAve_cut(repBoot), mseAve_cut(repBoot)] = computeMetrics(average(repBoot).betaCut(:,winERP),srcERP([1:4 13:18],winERP));
% [aucAveO_cut(repBoot), energyAveO_cut(repBoot), mseAveO_cut(repBoot)] = computeMetrics(overlap(repBoot).betaCut(:,winERP),srcERP([1:4 13:18],winERP));
% [aucAveNO_cut(repBoot), energyAveNO_cut(repBoot), mseAveNO_cut(repBoot)] = computeMetrics(noOverlap(repBoot).betaCut(:,winERP),srcERP([1:4 13:18],winERP));

end
end

figure;
subplot(1,2,1);hold on;
errorbar(1:4,[mean(mseAveO) mean(mseAveNO) mean(mseAveO_cut) mean(mseAveNO_cut)],[std(mseAveO) std(mseAveNO) std(mseAveO_cut) std(mseAveNO_cut)],'.')
bar(1:4,[mean(mseAveO) mean(mseAveNO) mean(mseAveO_cut) mean(mseAveNO_cut)])
ylabel('mse'); xticks(1:4); xticklabels({'overlap','no overlap','overlapCut','no overlap cut'})
subplot(1,2,2);hold on;
errorbar(1:4,[mean(energyAveO) mean(energyAveNO) mean(energyAveO_cut) mean(energyAveNO_cut)],[std(energyAveO) std(energyAveNO) std(energyAveO_cut) std(energyAveNO_cut)],'.')
bar(1:4,[mean(energyAveO) mean(energyAveNO) mean(energyAveO_cut) mean(energyAveNO_cut)])
ylabel('energy'); xticks(1:4); xticklabels({'overlap','no overlap','overlapCut','no overlap cut'})


%%% plot metrics
figure;
subplot(1,3,1);hold on;
errorbar(log(SNRlevel),mean(aucAve),std(aucAve),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveO),std(aucAveO),'LineWidth',2)
errorbar(log(SNRlevel),mean(aucAveNO),std(aucAveNO),'LineWidth',2)
xlabel('log(SNR)')
ylim([0 1])
ylabel('AUC')


clear u s v
[ u, s, v] = svd(sharedVal(5).avMapNoOverlap);
[ u1, s1, v1] = svd(sharedVal(3).avMapNoOverlap);
[ u2, s2, v2] = svd(sharedVal(6).avMapNoOverlap);
[ u3, s3, v3] = svd(sharedVal(6).avMapOverlap);
[ u4, s4, v4] = svd(sharedVal(10).avMapNoOverlap);
[ u5, s5, v5] = svd(sharedVal(10).avMapOverlap);
figure; 
subplot(3,2,1); plot( diag(s)./ sum(diag(s)) );title('fail avMapNoOverlap')
subplot(3,2,2); plot( diag(s1)./ sum(diag(s1)) );title('lowPerf avMapNoOverlap')
subplot(3,2,3); plot( diag(s2)./ sum(diag(s2)) );title('fail avMapNoOverlap')
subplot(3,2,4); plot( diag(s3)./ sum(diag(s3)) );title('fail avMapOverlap')
subplot(3,2,5); plot( diag(s4)./ sum(diag(s4)) );title('highPerf avMapNoOverlap')
subplot(3,2,6); plot( diag(s5)./ sum(diag(s5)) );title('highPerf avMapOverlap')
% figure; 
% subplot(3,2,1); plot( diag(s) );title('fail avMapNoOverlap')
% subplot(3,2,2); plot( diag(s1) );title('lowPerf avMapNoOverlap')
% subplot(3,2,3); plot( diag(s2) );title('fail avMapNoOverlap')
% subplot(3,2,4); plot( diag(s3) );title('fail avMapOverlap')
% subplot(3,2,5); plot( diag(s4) );title('highPerf avMapNoOverlap')
% subplot(3,2,6); plot( diag(s5) );title('highPerf avMapOverlap')


% condition number
testCond = diag(s);
testCond(1,1)/testCond(end,end)
cond(sharedVal(5).avMapNoOverlap)
cond(sharedVal(3).avMapNoOverlap)
cond(sharedVal(6).avMapNoOverlap)
cond(sharedVal(6).avMapOverlap)
cond(sharedVal(10).avMapNoOverlap)
cond(sharedVal(10).avMapOverlap)
% condition number is around 230

figure;plot(overlap(6).gcvError,'LineWidth',2)
hold on
plot(noOverlap(6).gcvError,'LineWidth',2)
plot(overlap(10).gcvError,'LineWidth',2)
plot(noOverlap(10).gcvError,'LineWidth',2)
legend('lowPerfSameSet [lambda=9977]','lowPerfDiffSet [lambda=1277]','highPerfSameSet [lambda=12597]','highPerfDiffSet [lambda=6552]')
title('gcvError')


figure;
for iSubj = 1:25
    subplot(5,5,iSubj);imagesc(sharedVal(6).roiMapOverlap(:,:,iSubj));
end
figure;
for iSubj = 1:25
    clear u s v
    [ u, s, v] = svd(sharedVal(10).roiMapOverlap(:,:,iSubj));
    subplot(5,5,iSubj);plot( diag(v)./trace(v) )
end

figure;
subplot(3,2,1);imagesc(sharedVal(5).avMapNoOverlap);title('fail avMapNoOverlap');
subplot(3,2,2);imagesc(sharedVal(3).avMapOverlap);title('lowPerf avMapNoOverlap')
subplot(3,2,3);imagesc(sharedVal(6).avMapNoOverlap);title('lowPerf avMapNoOverlap');
subplot(3,2,4);imagesc(sharedVal(6).avMapOverlap);title('lowPerf avMapOverlap');
subplot(3,2,5);imagesc(sharedVal(10).avMapNoOverlap);title('highPerf avMapNoOverlap');
subplot(3,2,6);imagesc(sharedVal(10).avMapOverlap);title('highPerf avMapOverlap');


figure;
for iSubj = 1:25
    subplot(5,5,iSubj);imagesc(sharedVal(7).fullFwd{iSubj});
end

figure;
subplot(2,2,1);imagesc(sharedVal(6).avMapNoOverlap);title('template low');colorbar;
subplot(2,2,2);imagesc(sharedVal(7).avMapNoOverlap);title('template low');colorbar;
subplot(2,2,3);imagesc(sharedVal(8).avMapNoOverlap);title('template higher');colorbar;
subplot(2,2,4);imagesc(sharedVal(5).avMapNoOverlap);title('template higher');colorbar;

testY = reshape(sharedVal(7).stackY,[size(Y_avg,2),size(Y_avg,1),size(Y_avg,3)]);
figure;
for iSubj = 1 :25
   subplot(5,5,iSubj);plotOnEgi(testY(:,iSubj,50))
end

figure;plotOnEgi(Y_avg(18,:,50))
figure;plotOnEgi(Y(18,:,50));colorbar



% cond number per participant?? 
cond(sharedVal(5).stackY)
cond(sharedVal(3).stackY)
cond(sharedVal(6).stackY)
cond(sharedVal(10).stackY)


[ uT, sT, vT] = csvd(sharedVal(5).avMapNoOverlap);

[bb, bmin, ll, gg, llmin] = minimum_norm(repmat(sharedVal(5).avMapNoOverlap,numSubs,1), sharedVal(5).stackY, 10);
[bbG, bminG, llG, ggG, llminG] = minimum_norm(repmat(sharedVal(5).avMapOverlap,numSubs,1), sharedVal(5).stackY, 10);
figure;imagesc(bb)
figure;imagesc(bbG)

figure;imagesc(noOverlap(5).beta)



%%% reduce ROI
cond(sharedVal(5).avMapOverlap(:,[1:4 13:18]))



%%%% test fitting
for iSub=1:50
    clear fwdMatrix roiInfo
    load([dataPath dirList(iSub).name])
    for rr=1:numROIs
        clear idxROI
        idxROI = find(strcmp(listROIs(rr),{roiInfo.name}));
        roiMapAll(:,rr,iSub) = sum(fwdMatrix(:,roiInfo(idxROI).meshIndices),2);
    end
end

for iSubj=1:50
    beta = avMap.\roiMapAll(:,:,iSubj);
    fittedSubj = beta.*avMap;
    leftAct(iSubj,:) = rms(fittedSubj - roiMapAll(:,:,iSubj));
end
figure;imagesc(leftAct)


%%%%%%%%%%%%%%%%
%%%% save matrices for Justin
templateNoOverlap = repmat(sharedVal(13).avMapNoOverlap,numSubs,1); 
templateOverlap = repmat(sharedVal(13).avMapOverlap,numSubs,1);
data = sharedVal(13).stackY;
figure;plot(data(:,80))
figure;plot(sharedVal(13).avMapOverlap(:,12))
minimum_norm(template, data, 10)
save('matrixMinNorm','data','templateOverlap','templateNoOverlap')


