% add results for template with best regularisation 
% in each of the simulation
addpath(genpath([pwd filesep 'subfunctions']))


%%%%% V1-MT simulation
clearvars; close all;
load ('simulOutput/simulV1MToutput.mat')
load('averageMap50Sum.mat') 

% get the variables from the data
nbBootstrap = size(simulERP,1);
sbjIncluded = zeros(1,size(simulERP,2));
for ss=1:size(simulERP,2)
    sbjIncluded(ss) = length(simulERP(1,ss,1).listSub);
end
noiseLevel = [simulERP(1,1,:).noise];

% do the loops
for bb=1:nbBootstrap
    for nbSbj = 1: length(sbjIncluded)
        fprintf('N%d bootstrap %d\n',sbjIncluded(nbSbj),bb)
        for noise = 1:length(noiseLevel)
            simulData = simulERP(bb,nbSbj,noise).data;
            sourceData = simulERP(bb,nbSbj,noise).srcERP;
            [beta, betaCurv, betaBest, lambda, lambdaCurv, lambdaBest, ...
                lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(simulData,1)),sourceData);
%             figure;imagesc(betaBest);figure;imagesc(betaCurv);figure;imagesc(beta)
            simulERP(bb,nbSbj,noise).beta(8,:,:) = betaBest;
        end
    end
end

save('simulOutput/simulV1MToutput.mat','simulERP','-v7.3')


%%%%% V1-MT simulation

clearvars; close all;
load ('simulOutput/simulV2V4output.mat')
load('averageMap50Sum.mat') 

% get the variables from the data
nbBootstrap = size(simulERP,1);
sbjIncluded = zeros(1,size(simulERP,2));
for ss=1:size(simulERP,2)
    sbjIncluded(ss) = length(simulERP(1,ss,1).listSub);
end
noiseLevel = [simulERP(1,1,:).noise];

% do the loops
for bb=1:nbBootstrap
    for nbSbj = 1: length(sbjIncluded)
        fprintf('N%d bootstrap %d\n',sbjIncluded(nbSbj),bb)
        for noise = 1:length(noiseLevel)
            simulData = simulERP(bb,nbSbj,noise).data;
            sourceData = simulERP(bb,nbSbj,noise).srcERP;
            [beta, betaCurv, betaBest, lambda, lambdaCurv, lambdaBest, ...
                lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(simulData,1)),sourceData);
%             figure;imagesc(betaBest);figure;imagesc(betaCurv);figure;imagesc(beta)
            simulERP(bb,nbSbj,noise).beta(8,:,:) = betaBest;
        end
    end
end

save('simulOutput/simulV2V4output.mat','simulERP','-v7.3')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% compare systems V1MT
clearvars; close all;
load ('simulOutput/simulSysV1MT.mat')
dirModel = '/Users/marleneponcet/Documents/Git/templateMinNorm/createTemplate/templates/';

% get the variables from the data
nbBootstrap = size(simulSys,1);
noiseLevel = [simulSys(1,:,1).noise];
nbSystem = size(simulSys,3);

% do the loops
for sys = 1:nbSystem
    clear avMap
    switch sys
        case 1
            load([dirModel 'averageMapEGI32.mat']);
        case 2
            load([dirModel 'averageMapEGI64.mat']);
        case 3
            load([dirModel 'averageMapEGI128.mat']);
        case 4
            load([dirModel 'averageMapEGI256.mat']);
    end
    for bb=1:nbBootstrap
        fprintf('model%d bootstrap %d\n',sys,bb)
        for noise = 1:length(noiseLevel)
            simulData = simulSys(bb,noise,sys).data;
            sourceData = simulSys(bb,noise,sys).srcERP;
            [beta, betaCurv, betaBest, lambda, lambdaCurv, lambdaBest, ...
                lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(simulData,1)),sourceData);
            %             figure;imagesc(betaBest);figure;imagesc(betaCurv);figure;imagesc(beta)
            simulSys(bb,noise,sys).beta(5,:,:) = betaBest;
        end
    end
end

save('simulOutput/simulSysV1MT.mat','simulSys','-v7.3')



%%%%% compare systems V2V4
clearvars; close all;
load ('simulOutput/simulSysV2V4.mat')
dirModel = '/Users/marleneponcet/Documents/Git/templateMinNorm/createTemplate/templates/';

% get the variables from the data
nbBootstrap = size(simulSys,1);
noiseLevel = [simulSys(1,:,1).noise];
nbSystem = size(simulSys,3);

% do the loops
for sys = 1:nbSystem
    clear avMap
    switch sys
        case 1
            load([dirModel 'averageMapEGI32.mat']);
        case 2
            load([dirModel 'averageMapEGI64.mat']);
        case 3
            load([dirModel 'averageMapEGI128.mat']);
        case 4
            load([dirModel 'averageMapEGI256.mat']);
    end
    for bb=1:nbBootstrap
        fprintf('model%d bootstrap %d\n',sys,bb)
        for noise = 1:length(noiseLevel)
            simulData = simulSys(bb,noise,sys).data;
            sourceData = simulSys(bb,noise,sys).srcERP;
            [beta, betaCurv, betaBest, lambda, lambdaCurv, lambdaBest, ...
                lambdaGridMinNorm] = minNorm_lcurve_bestRegul(avMap, squeeze(mean(simulData,1)),sourceData);
            %             figure;imagesc(betaBest);figure;imagesc(betaCurv);figure;imagesc(beta)
            simulSys(bb,noise,sys).beta(5,:,:) = betaBest;
        end
    end
end

save('simulOutput/simulSysV2V4.mat','simulSys','-v7.3')

