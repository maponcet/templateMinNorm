clearvars;close all;
% compare ERP for different noise and nb of sbj with the same source amplitude
% be aware that for N=1 it is NOT the same sbj that is plotted for
% different SNR!

addpath(genpath([pwd filesep 'subfunctions']))
dataPath = '/Users/marleneponcet/Documents/data/skeriDATA/forwardAllEGI/';
dirList = dir([dataPath 'forward*']);
load('averageMap50Sum.mat') % load average map of ROIs (128 elec x 18 ROIs)
numROIs = length(listROIs);

% some parameters
SNRlevel = [0.1 1 10 200 10000]; % noise level 10% = if the signal is 10 then the noise is 10*10, for 0.5 S=50/100
nLambdaRidge = 50; % for calculating minimum_norm, reg constant, hyper param in min norm
numCols = 5; % For reducing dimensionality of data: use first X columns of v ([~, ~, v] = svd(Y);) as time basis (old code = 2, new = 5)
% set 2 vectors of the left and right sources in the same order
sourceL = {'V1-L','MT-L'};
sourceR = {'V1-R','MT-R'};
% simulated signal
activeROIs = [sourceL,sourceR]; % left sources then right sources to make it easier for copying the same signal over the 2 hemispheres
% find the ROI index corresponding to the activeROIs
ac_sources = cell2mat(arrayfun(@(x) cellfind(listROIs,activeROIs{x}),1:length(activeROIs),'uni',false));

% nbSbjToInclude =[1 2 5 10 20 30 40 50];
numSubs = 50;
snrRawTime = zeros(5,length(SNRlevel));


        
for test=1:5   
    
%% Simulate sources (sourceERP)
% amplitude (1 to 10) and time function is different for each
% source but the same for all sbj for a given bootstrap
[srcAmp, srcSSVEP, srcERP,winERP] = createSourceROI(numROIs,ac_sources(1:length(ac_sources)/2),ac_sources((length(ac_sources)/2+1):end));
% ERP baseline timewindow
timeBase = setdiff(1:size(srcERP,2),winERP);

% clf;
f1=figure; set(gcf,'position',[100,100,1200,700]);hold on;
% f2=figure; set(gcf,'position',[100,100,1200,700]);hold on;


% list of random sbj with replacement
listSub = randi(length(dirList),numSubs,1);

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
        
        
for level=1:length(SNRlevel)
        noiseLevel = SNRlevel(level);           
          
        %% Simulate scalp activity (Y)
        % use the generated sources to simulate scalp activity for each sbj 
        % (using individual fwd model)
        Y = zeros(numSubs,size(fullFwd{1},1),length(srcERP));
        Ylo = Y; Y_noise = Y; 
        
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
                
        %% Use average reference for centering Y 
        for iSub=1:numSubs
            Y(iSub,:,:) = bsxfun(@minus,squeeze(Y(iSub,:,:)), mean(squeeze(Y(iSub,:,:))));
%             Y2(iSub,:,:) = bsxfun(@minus,squeeze(Y2(iSub,:,:)), mean(squeeze(Y2(iSub,:,:))));
        end
 
        % check SNR
        ySNR_time=zeros(1,numSubs);
        for iSub=1:numSubs
            ySNR_time(iSub) = (rms(rms(Y(iSub,:,winERP)))/rms(rms(Y(iSub,:,timeBase)))) ^2 -1;
        end       
        snrRawTime(test,level) = mean(ySNR_time); 
        
%         figure(f1);
%         subplot(2,3,level);hold on;
%         plot(squeeze(mean(Y(1:2,1,:),1))); % plot average 2 sbj
%         plot(squeeze(mean(Y(1:8,1,:),1))); % plot average 8 sbj
%         plot(squeeze(mean(Y(1:20,1,:),1))); % plot average 20 sbj
%         plot(squeeze(mean(Y(1:50,1,:),1))); % plot average 50 sbj        
%         legend('N=2','N=8','N=20','N=50')
%         title(['SNR' num2str(SNRlevel(level))])
%         figure(f2); 
%         subplot(2,3,level);hold on;
%         for iSub = 1:4
%             plot(squeeze(Y(iSub,1,:))); % plot S1
%         end
%         title(['SNR' num2str(SNRlevel(level))])
%         legend('S1','S2','S3','S4')
            
        figure(f1);
        subplot(2,5,level);hold on;
        plot(squeeze(Y(15,18,:))); 
        plot(squeeze(Y(5,18,:))); 
        plot(squeeze(Y(10,18,:))); 
        legend('S2','S5','S10')
        title(['SNR ' num2str(SNRlevel(level))])
        subplot(2,5,level+5);hold on;
        plot(squeeze(mean(Y(1:20,18,:),1)));
        legend('average 20')
        
end % SNR
% saveas(f1,['figures' filesep 'checkSNRt' num2str(test)],'png')
% saveas(f2,['figures' filesep 'checkSNRindT' num2str(test)],'png')
saveas(f1,['figures' filesep 'checkSNRmixT' num2str(test)],'png')
end

