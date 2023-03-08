%% %%%%%%%%%%%%%%%%%% PRECEDING TRIALS (Gambler's fallacy effect) %%%%%%%%%%%%%%%%%%%%%%

clear; close all;clc
mainDir = 'D:\presentiment_eeg';
dataDir = fullfile(mainDir, 'data_raw');
codeDir = fullfile(mainDir, 'paa_eeg');
load(fullfile(codeDir, 'sInfo_old2.mat'));

% number of preceding trials to examine
ntrials = 2;  

% Gather 3 markers preceding each event
pleasant = [];
neutral = [];
unpleasant = [];
emotional = [];
progressbar('Getting event info')
for iSub = 1:length(sInfo)
    fprintf('Subject %2.2d \n', iSub)
    
    if sInfo(iSub).group == 1
        filepath = fullfile(dataDir, 'group-01');
    else
        filepath = fullfile(dataDir, 'group-02');
    end
    EEG = pop_loadcnt(fullfile(filepath, [sInfo(iSub).filename '.cnt']), ...
        'dataformat','auto','keystroke','on','memmapfile','');

    events = str2double({EEG.event.type});
    events(isnan(events)) = [];
    events(events == 255 | events == 239 | events == 223 | events == 207) = [];
    summary(categorical(events))
    error('check that events are not already numbers before checking randomness!')
    [~,p_random(iSub)] = runstest(events);
    if p_random(iSub) < 0.05, warning('This trial sequence is no random'); end

%     pleasant = [];
%     neutral = [];
%     unpleasant = [];

    for iEv = 4:length(events)
        % Pleasant condition
        if events(iEv) == 2
            pleasant = [ pleasant events(iEv-ntrials:iEv-1) ];
%             count1 = count1 + 1;

            % track Pleasant trials preceded by 3 neutral ones
            if sum( events(iEv-ntrials:iEv-1) == 4 ) == 3
                expectPleasant(iEv) = 1;
            else
                expectPleasant(iEv) = 0;
            end

        % Neutral condition
        elseif events(iEv) == 4
            neutral = [ neutral events(iEv-ntrials:iEv-1) ];

            % track neutral trials with 3 previous emotional trials
            if sum( events(iEv-ntrials:iEv-1) == 2 | events(iEv-ntrials:iEv-1) == 8 ) == 3
                expectNeutral(iEv) = 1;
            else
                expectNeutral(iEv) = 0;
            end

%             count2 = count2 + 1;
        
        % Unpleasant condition
        elseif events(iEv) == 8
            unpleasant = [ unpleasant events(iEv-ntrials:iEv-1) ];
%             count3 = count3 + 1;

            % track Unplesant trials preceded by 3 neutral ones
            if sum( events(iEv-ntrials:iEv-1) == 4 ) == 3
                expectUnpleasant(iEv) = 1;
            else
                expectUnpleasant(iEv) = 0;
            end

        % Emotional condition
        elseif events(iEv) == 2 || events(iEv) == 8

            % track emotional trials preceded by 3 neutral trials
            emotional =  [ emotional events(iEv-ntrials:iEv-1) ];
            if sum( events(iEv-ntrials:iEv-1) == 4 ) == 3
                expectEmotional(iEv) = 1;
            else
                expectEmotional(iEv) = 0;
            end
        end

        progressbar([], iEv/length(events));
    end
    
%     subplot(3,1,1); 
%     histogram(pleasant); hold on;
%     subplot(3,1,2); 
%     histogram(neutral); hold on;
%     subplot(3,1,3); 
%     histogram(unpleasant); hold on;

%     Pleasant(:,iSub) = { pleasant };
%     Neutral(:,iSub) = { neutral };
%     Unpleasant(:,iSub) = { unpleasant };

    progressbar(iSub/length(sInfo));
end

% subplot(3,1,1); title('3 events preceding Pleasant trials')
% subplot(3,1,2); title('3 events preceding Neutral trials')
% subplot(3,1,3); title('3 events preceding Unpleasant trials')

% % Remove bad events
% pleasant(isnan(pleasant)) = [];
% neutral(isnan(neutral)) = [];
% unpleasant(isnan(unpleasant)) = [];
% pleasant(pleasant == 255 | pleasant == 239 | pleasant == 223 | pleasant == 207) = [];
% neutral(neutral == 255 | neutral == 239 | neutral == 223 | neutral == 207) = [];
% unpleasant(unpleasant == 255 | unpleasant == 239 | unpleasant == 223 | unpleasant == 207) = [];

subplot(3,1,1);
histogram(pleasant); hold on;
subplot(3,1,2);
histogram(neutral); hold on;
subplot(3,1,3);
histogram(unpleasant); hold on;

% summary(categorical(pleasant))
% summary(categorical(neutral))
% summary(categorical(unpleasant))

% figure;
% subplot(1,3,1)
% boxplot(pleasant); 
% subplot(1,3,2)
% boxplot(neutral)
% subplot(1,3,3)
% boxplot(unpleasant)

% Run stats
% for iSub = 1:length(sInfo)
%     x1 = Neutral{iSub};
%     x2 = Pleasant{iSub};
%     h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
%     if h1==0 && h2==0 %&& h3==0
%         cprintf('blue','Normal distribution --> two-sample t-test. \n');
%         [~,df(iSub),CI(iSub,:),~,~,~,p(iSub)]= limo_ttest(2,x1,x2,0.05);
%     else
%          cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
%         %  [p(iFrame),h,stats] = signrank(y,X)
%         [t(iSub),diff(iSub),CI(iSub,:),~,~,df(iSub),p(iSub)] = yuen(x1,x2,20);
%     end
% end
% [h, ~, ~, adj_p] = fdr_bh(p,0.05,'pdep','yes');


%% Simulations: count number of trials preceded by 3 consecutive
nboot = 10000;
for iboot = 1:nboot
    
    for iTrial = 1:size(data,1)
        tmpidx = ones(1,3);  % index of options 
        tmpidx(guesser(iTrial)) = 0; % remove guesser from options as it can't guess itself
        caller_H0(iTrial,iboot) = randsample(find(tmpidx),1);   % random caller
    end
end

 
save(fullfile(outDir, '3prev_trials.mat'), 'results')
save(fullfile(outDir, '3prev_trials_H0.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
sum(mask)


%% %%%%%%%%%%%%%%%%%% TIME-ON-TASK EFFECT %%%%%%%%%%%%%%%%%%%

clear; close all;clc
dataDir = 'D:\presentiment_eeg\data_clean';
outDir = 'D:\presentiment_eeg\order_effects';
cd(dataDir)
eeglab; close;

[STUDY, ALLEEG] = pop_loadstudy('filename','presentiment_eeg.study','filepath',dataDir);
disp('done')

%% Compare conditions under H0 across all trials from raw ERP data

pleasant = [];
neutral = [];
unpleasant = [];

for iSub = 1:nSub
    fprintf('Subject %2.2d \n', iSub)
    tmp = load('-mat', fullfile(STUDY.datasetinfo(iSub).filepath, sprintf('sub-%2.2d.daterp',iSub)));
        
    times = tmp.times;
    xAxis = times > -1000 & times < 10;
    
    conds = [tmp.trialinfo.type]';
    idxCond1 = strcmp(string(conds), '2');
    idxCond2 = strcmp(string(conds), '4');
    idxCond3 = strcmp(string(conds), '8');
    
    tmp = rmfield(tmp, { 'labels' 'times' 'datatype' 'parameters' 'datafiles' 'datatrials' 'trialinfo'});
    tmp = struct2table(tmp);
    
    for iChan = 1:64
        erpData = tmp{:,iChan};   
        pleasant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond1),20,2);
        neutral(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond2),20,2);
        unpleasant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond3),20,2);
    end
end

% Get channel neighbors for cluster correction for multiple comparison
chanlocs = ALLEEG(1).chanlocs;
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs,0);

% Pleasant-Neutral under H0 (all trials)
[results, results_H0] = compute_randomeffect(pleasant,neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'all_pairwise_comparisons.mat'), 'results')
save(fullfile(outDir, 'all_pairwise_comparisons_H0.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('All trials (uncorrected)')
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('All trials (cluster-corrected)')

% % Unpleasant-Neutral under H0 (all trials)
% [results, results_H0] = compute_randomeffect(unpleasant,neutral,1000,'trimmed mean');
% tvals = results(:,:,1);
% pvals = results(:,:,2);
% tvals_H0 = squeeze(results_H0(:,:,1,:));
% pvals_H0 = squeeze(results_H0(:,:,2,:));
% save(fullfile(outDir, 'all_pairwise_comparisons_unpleasant.mat'), 'results')
% save(fullfile(outDir, 'all_pairwise_comparisons_H0_unpleasant.mat'), 'results_H0', '-v7.3')
% mcctype = 0;
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
% plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
% title('All trials (uncorrected)')
% mcctype = 2;
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
% plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
% title('All trials (cluster-corrected)')

%% Compare across 3 time blocks from raw ERP data

% Calculate minimum block size
tmp = { STUDY.datasetinfo.trialinfo }';
nSub = length(tmp);
for i = 1:nSub
    trialLength(i,:) = length(STUDY.datasetinfo(i).trialinfo);
end
minTrial = min(trialLength);
blockSize = floor(minTrial / 3);
clear tmp

% chan = 64;   % 64 = O2

% [STUDY, erpdata, erptimes] = std_erpplot(STUDY,ALLEEG,'channels',{EEG.chanlocs(chan).labels},'timerange',[-1000 10],'noplot','on');

for iSub = 1:nSub
    fprintf('Subject %2.2d \n', iSub)
    tmp = load('-mat', fullfile(STUDY.datasetinfo(iSub).filepath, sprintf('sub-%2.2d.daterp',iSub)));
        
    times = tmp.times;
    xAxis = times > -1000 & times < 10;
    
    conds = [tmp.trialinfo.type]';
    
    tmp = rmfield(tmp, { 'labels' 'times' 'datatype' 'parameters' 'datafiles' 'datatrials' 'trialinfo'});
    tmp = struct2table(tmp);
    
    % Trial block index and count number of trials
    % block 1
    idx = 1:blockSize;
    idxCond1 = strcmp(string(conds(idx,:)), '2');
    idxCond2 = strcmp(string(conds(idx,:)), '4');
    idxCond3 = strcmp(string(conds(idx,:)), '8');
    block1_nTrials(iSub,1) = sum(idxCond1);
    block1_nTrials(iSub,2) = sum(idxCond2);
    block1_nTrials(iSub,3) = sum(idxCond3);
    % fprintf('Pleasant:    %d trials \n', sum(idxCond1))
    % fprintf('Neutral:     %d trials \n', sum(idxCond2))
    % fprintf('Unpleasant:  %d trials \n', sum(idxCond3))
    for iChan = 1:64
        erpData = tmp{:,iChan};   
        block1_pleasant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond1),10,2);
        block1_neutral(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond2),10,2);
        block1_unplesant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond3),10,2);
    end

    % block 2
    idx = blockSize+1:blockSize*2;
    idxCond1 = strcmp(string(conds(idx,:)), '2');
    idxCond2 = strcmp(string(conds(idx,:)), '4');
    idxCond3 = strcmp(string(conds(idx,:)), '8');
    block2_nTrials(iSub,1) = sum(idxCond1);
    block2_nTrials(iSub,2) = sum(idxCond2);
    block2_nTrials(iSub,3) = sum(idxCond3);
    for iChan = 1:64
        erpData = tmp{:,iChan};   
        block2_pleasant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond1),10,2);
        block2_neutral(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond2),10,2);
        block2_unplesant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond3),10,2);
    end

    % block 3
    idx = blockSize*2+1:trialLength(iSub);
    idxCond1 = strcmp(string(conds(idx,:)), '2');
    idxCond2 = strcmp(string(conds(idx,:)), '4');
    idxCond3 = strcmp(string(conds(idx,:)), '8');
    block3_nTrials(iSub,1) = sum(idxCond1);
    block3_nTrials(iSub,2) = sum(idxCond2);
    block3_nTrials(iSub,3) = sum(idxCond3);
    for iChan = 1:64
        erpData = tmp{:,iChan};   
        block3_pleasant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond1),10,2);
        block3_neutral(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond2),10,2);
        block3_unplesant(iChan,:,iSub) = trimmean(erpData(xAxis,idxCond3),10,2);
    end
end

% Compare number of trials per condition per block
% figure('color','w');
% subplot(3,1,1)
% histogram(block1_nTrials(:,1)); hold on; histogram(block1_nTrials(:,2)); hold on; histogram(block1_nTrials(:,3));
% legend('Pleasant', 'Neutral', 'Unpleasant')
% subplot(3,1,2)
% histogram(block2_nTrials(:,1)); hold on; histogram(block2_nTrials(:,2)); hold on; histogram(block2_nTrials(:,3));
% legend('Pleasant', 'Neutral', 'Unpleasant')
% subplot(3,1,3)
% histogram(block3_nTrials(:,1)); hold on; histogram(block3_nTrials(:,2)); hold on; histogram(block3_nTrials(:,3));
% legend('Pleasant', 'Neutral', 'Unpleasant')

% Get channel neighbors for cluster correction for multiple comparison
chanlocs = ALLEEG(1).chanlocs;
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs,0);

disp('done') 

% Compare conditions under H0 for Block 1
[results, results_H0] = compute_randomeffect(block1_pleasant,block1_neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'block1_pairwise_comparisons.mat'), 'results')
save(fullfile(outDir, 'block1_pairwise_comparisons_H0.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('Block 1')
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('Block 1 (cluster-corrected)')
gong

% Compare conditions under H0 for Block 2
[results, results_H0] = compute_randomeffect(block2_pleasant,block2_neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'block2_pairwise_comparisons.mat'), 'results')
save(fullfile(outDir,'block2_pairwise_comparisons_H0.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('Block 2')
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('Block 2 (cluster-corrected)')
gong

% Compare conditions under H0 for Block 3
[results, results_H0] = compute_randomeffect(block3_pleasant,block3_neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'block3_pairwise_comparisons.mat'), 'results')
save(fullfile(outDir,'block3_pairwise_comparisons_H0.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('Block 3')
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', times(xAxis), tvals, mask, pcorr, mcctype, chanlocs)
title('Block 3 (cluster-corrected)')
gong






%% Compare across all trials on LIMO files

cd 'D:\presentiment_eeg\data_clean\derivatives\LIMO_presentiment_eeg\pleasant_neutral'
tmp = load('LIMO.mat');
gpLIMO = tmp.LIMO;
% times = gpLIMO.data.timevect;
% peakTime = find(times == -148);
tmp = { STUDY.datasetinfo.trialinfo }';
nSub = length(tmp); clear tmp

pleasant = nan(64,330,78);
neutral = nan(64,330,78);
unpleasant = nan(64,330,78);

% figure('color','w'); 
progressbar('Gathering mean data for each subject and condition')
for iSub = 1:nSub
    fprintf('Subject %2.2d \n', iSub)
    subPath = gpLIMO.data.data{1, 1}{iSub};
    subPath = subPath(1:end-10);
    load(fullfile(subPath, 'LIMO.mat'))     % Structure containing all data and model information
    conds = LIMO.data.Cat; 
    xAxis = LIMO.data.timevect';  
%     load(fullfile(subPath, 'Betas.mat'))    % Beta parameter estimates
    load(fullfile(subPath, 'Yr.mat'))       % Single trial data reorganized to fit X (grouped by condition)
%     load(fullfile(subPath, 'Yhat.mat'))     % Predicted data 
%     load(fullfile(subPath, 'Res.mat'))      % Residuals (non modelled) data 
%     load(fullfile(subPath, 'Condition_effect_1.mat')) % Factor effect (Categorical designs)
%     fvals = Condition_effect(64,:,1);
%     pvals = Condition_effect(64,:,2);
%     sigLat = find(pvals < 0.05);
%     xAxis(sigLat)
%     
%     hold on;
%     for iTrial = 1:size(Yr,3)
%         if conditions(iTrial) == 1
%             plot(xAxis,squeeze(Yr(64,:,iTrial)),'g'); 
%             plot(xAxis,squeeze(Yhat(64,:,iTrial)),'g'); 
%             plot(xAxis,squeeze(Res(64,:,:)),'g'); 
%         elseif conditions(iTrial) == 2
%             plot(xAxis,squeeze(Yr(64,:,iTrial)),'k'); 
%             plot(xAxis,squeeze(Yhat(64,:,iTrial)),'k'); 
%             plot(xAxis,squeeze(Res(64,:,:)),'k'); 
%         elseif conditions(iTrial) == 3
%             plot(xAxis,squeeze(Yr(64,:,iTrial)),'r'); 
%             plot(xAxis,squeeze(Yhat(64,:,iTrial)),'r'); 
%             plot(xAxis,squeeze(Res(64,:,:)),'r'); 
%         end
%     end
    
    for iChan = 1:64
%         erpData = squeeze(Res(iChan,:,:));   % on Res data 
        erpData = squeeze(Yr(iChan,:,:));   % on Yr data 
%         erpData = squeeze(Yhat(iChan,:,:));   % on Yhat data 
        pleasant(iChan,:,iSub) = trimmean(erpData(:,conds == 1),20,2);
        neutral(iChan,:,iSub) = trimmean(erpData(:,conds == 2),20,2);
        unpleasant(iChan,:,iSub) = trimmean(erpData(:,conds == 3),20,2);
    end

    progressbar(iSub / nSub)
end
% legend('Pleasant', 'Neutral', 'Unpleasant')

% Get channel neighbors for cluster correction for multiple comparison
chanlocs = ALLEEG(1).chanlocs;
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs,0);

% Pleasant-Neutral under H0 (all trials)
[results, results_H0] = compute_randomeffect(pleasant,neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'all_pairwise_comparisons_LIMO.mat'), 'results')
save(fullfile(outDir, 'all_pairwise_comparisons_H0_LIMO.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
% title('All trials (uncorrected)')
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
title('All trials (cluster-corrected)')

% % Unpleasant-Neutral under H0 (all trials)
% [results, results_H0] = compute_randomeffect(unpleasant,neutral,1000,'trimmed mean');
% tvals = results(:,:,1);
% pvals = results(:,:,2);
% tvals_H0 = squeeze(results_H0(:,:,1,:));
% pvals_H0 = squeeze(results_H0(:,:,2,:));
% save(fullfile(outDir, 'all_pairwise_comparisons_unpleasant_LIMO.mat'), 'results')
% save(fullfile(outDir, 'all_pairwise_comparisons_H0_unpleasant_LIMO.mat'), 'results_H0', '-v7.3')
% mcctype = 0;
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
% plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
% % title('All trials (uncorrected)')
% mcctype = 2;
% [mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
% plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
% title('All trials (cluster-corrected)')

gong

%% Compare across all trials on Betas to replicate effect

pleasant = nan(64,330,1);
neutral = nan(64,330,1);
unpleasant = nan(64,330,1);

for iSub = 1:nSub
    fprintf('Subject %2.2d \n', iSub)
    subPath = gpLIMO.data.data{1, 1}{iSub};
    subPath = subPath(1:end-10);
    load(fullfile(subPath, 'LIMO.mat'))     % Structure containing all data and model information
    conds = LIMO.data.Cat; 
    xAxis = LIMO.data.timevect';  
    load(fullfile(subPath, 'Betas.mat'))    % Beta parameter estimates
    
    for iChan = 1:64
        pleasant(iChan,:,iSub) = squeeze(Betas(iChan,:,1));
        neutral(iChan,:,iSub) = squeeze(Betas(iChan,:,2));
        unpleasant(iChan,:,iSub) = squeeze(Betas(iChan,:,3));
    end
end

% Get channel neighbors for cluster correction for multiple comparison
chanlocs = ALLEEG(1).chanlocs;
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs,0);

% Pleasant-Neutral under H0 (all trials)
[results, results_H0] = compute_randomeffect(pleasant,neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'all_pairwise_comparisons_LIMO_Betas.mat'), 'results')
save(fullfile(outDir, 'all_pairwise_comparisons_H0_LIMO_Betas.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
% title('All trials (uncorrected)')
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
title('All trials (cluster-corrected)')

gong

%% Compare across 3 blocks from LIMO files

pleasant = nan(64,330,78);
neutral = nan(64,330,78);
unpleasant = nan(64,330,78);

progressbar('Gathering data for each subject and condition')
for iSub = 1:nSub
    fprintf('Subject %2.2d \n', iSub)
    subPath = gpLIMO.data.data{1, 1}{iSub};
    subPath = subPath(1:end-10);
    load(fullfile(subPath, 'LIMO.mat'))     % Structure containing all data and model information
    conds = LIMO.data.Cat; 
    xAxis = LIMO.data.timevect';  
%     load(fullfile(subPath, 'Betas.mat'))    % Beta parameter estimates
%     load(fullfile(subPath, 'Yr.mat'))       % Single trial data reorganized to fit X (grouped by condition)
%     load(fullfile(subPath, 'Yhat.mat'))     % Predicted data 
    load(fullfile(subPath, 'Res.mat'))      % Residuals (non modelled) data 
%     load(fullfile(subPath, 'Condition_effect_1.mat')) % Factor effect (Categorical designs)
%     fvals = Condition_effect(64,:,1);
%     pvals = Condition_effect(64,:,2);
%     sigLat = find(pvals < 0.05);
%     xAxis(sigLat)
%     
%     figure('color','w'); hold on;
%     for iTrial = 1:size(Yr,3)
%         if conditions(iTrial) == 1
%             plot(xAxis,squeeze(Yr(64,:,iTrial)),'g'); 
%             plot(xAxis,squeeze(Yhat(64,:,iTrial)),'g'); 
%             plot(xAxis,squeeze(Res(64,:,:)),'g'); 
%         elseif conditions(iTrial) == 2
%             plot(xAxis,squeeze(Yr(64,:,iTrial)),'k'); 
%             plot(xAxis,squeeze(Yhat(64,:,iTrial)),'k'); 
%             plot(xAxis,squeeze(Res(64,:,:)),'k'); 
%         elseif conditions(iTrial) == 3
%             plot(xAxis,squeeze(Yr(64,:,iTrial)),'r'); 
%             plot(xAxis,squeeze(Yhat(64,:,iTrial)),'r'); 
%             plot(xAxis,squeeze(Res(64,:,:)),'r'); 
%         end
%     end
    
    for iChan = 1:64
        erpData = squeeze(Res(iChan,:,:));   % on Res data 
        pleasant(iChan,:,iSub) = trimmean(erpData(:,conds == 1),20,2);
        neutral(iChan,:,iSub) = trimmean(erpData(:,conds == 2),20,2);
        unpleasant(iChan,:,iSub) = trimmean(erpData(:,conds == 3),20,2);
    end

    progressbar(iSub / nSub)
end

% Get channel neighbors for cluster correction for multiple comparison
chanlocs = ALLEEG(1).chanlocs;
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs,0);

% Pleasant-Neutral under H0 (all trials)
[results, results_H0] = compute_randomeffect(pleasant,neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'all_pairwise_comparisons_LIMO.mat'), 'results')
save(fullfile(outDir, 'all_pairwise_comparisons_H0_LIMO.mat'), 'results_H0', '-v7.3')
mcctype = 0;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
% title('All trials (uncorrected)')
mcctype = 2;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
plotresults('time', xAxis, tvals, mask, pcorr, mcctype, chanlocs)
title('All trials (cluster-corrected)')





%% Check WLS weights

cd 'D:\presentiment_eeg\data_clean\derivatives\LIMO_presentiment_eeg\pleasant_neutral'
tmp = load('LIMO.mat');
gpLIMO = tmp.LIMO;

% Create list of LIMO.mat
subPath = {};
for iSub = 1:nSub
    subPath(iSub) = { [tmp.LIMO.data.data_dir{iSub} '\LIMO.mat'] };
end
chanlocpath = 'D:\presentiment_eeg\data_clean\derivatives\limo_gp_level_chanlocs.mat';
limo_CheckWeight(subPath, chanlocpath,'CheckBias','on', ...
    'TestDifference','on', 'SingleSubjectAnalysis','on','PlotRank','on')


gong



%% Stats comparing across conditions for each block (not under H0)

% ERP BLOCK 1
figure('color','w');
for iFrame = 1:size(block1_pleasant,1)
    x1 = block1_pleasant(iFrame,:);
    x2 = block1_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        % cprintf('blue','Normal distribution --> two-sample t-test. \n');
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(3,1,1)
plotDiff(times(xAxis),block1_pleasant,block1_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Block 1')

% ERP BLOCK 2
for iFrame = 1:size(block2_pleasant,1)
    x1 = block2_pleasant(iFrame,:);
    x2 = block2_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        % cprintf('blue','Normal distribution --> two-sample t-test. \n');
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(3,1,2)
plotDiff(times(xAxis),block2_pleasant,block2_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Block 2')

% ERP BLOCK 3
for iFrame = 1:size(block3_pleasant,1)
    x1 = block3_pleasant(iFrame,:);
    x2 = block3_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        % cprintf('blue','Normal distribution --> two-sample t-test. \n');
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(3,1,3)
plotDiff(times(xAxis),block3_pleasant,block3_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Block 3')

%% Try stats on all trials to see if there is an effect
pleasant = [];
neutral = [];
unpleasant = [];

for iSub = 1:nSub
    fprintf('Subject %2.2d \n', iSub)
    tmp = load('-mat', fullfile(STUDY.datasetinfo(iSub).filepath, sprintf('sub-%2.2d.daterp',iSub)));

    erpData = tmp.chan64;   % O2
    times = tmp.times;
    xAxis = times > -1000 & times < 1000;
    conds = [tmp.trialinfo.type]';

    idxCond1 = strcmp(string(conds(idx,:)), '2');
    idxCond2 = strcmp(string(conds(idx,:)), '4');
    idxCond3 = strcmp(string(conds(idx,:)), '8');

    pleasant(:,iSub) = trimmean(erpData(xAxis,idxCond1),20,2);
    neutral(:,iSub) = trimmean(erpData(xAxis,idxCond2),20,2);
    unpleasant(:,iSub) = trimmean(erpData(xAxis,idxCond3),20,2);

end

% Pleasant-neutral
figure('color','w');
subplot(2,1,1)
for iFrame = 1:size(pleasant,1)
    x1 = pleasant(iFrame,:);
    x2 = neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2);
    if h1==0 && h2==0
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
h = fdr_bh(p,0.05,'pdep','no');
plotDiff(times(xAxis),pleasant,neutral,'trimmed mean', h,'Pleasant','Neutral')
title('All trials (channel O2)')

% Pleasant-neutral
subplot(2,1,2)
for iFrame = 1:size(pleasant,1)
    x1 = unpleasant(iFrame,:);
    x2 = neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2);
    if h1==0 && h2==0
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
h = fdr_bh(p,0.05,'pdep','no');
plotDiff(times(xAxis),unpleasant,neutral,'trimmed mean', h,'Unpleasant','Neutral')
title('All trials (channel O2)')

%% Comparing same condition across blocks

% Pleasant BLOCK 1-2
figure('color','w');
for iFrame = 1:size(block1_pleasant,1)
    x1 = block1_pleasant(iFrame,:);
    x2 = block2_pleasant(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,1)
plotDiff(times(xAxis),block1_pleasant,block2_pleasant,'trimmed mean', h,'Pleasant','Neutral')
title('Pleasant Block 1-2')

% Pleasant BLOCK 2-3
for iFrame = 1:size(block1_pleasant,1)
    x1 = block2_pleasant(iFrame,:);
    x2 = block3_pleasant(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,2)
plotDiff(times(xAxis),block2_pleasant,block3_pleasant,'trimmed mean', h,'Pleasant','Neutral')
title('Pleasant Block 2-3')

% Neutral BLOCK 1-2
for iFrame = 1:size(block1_pleasant,1)
    x1 = block1_neutral(iFrame,:);
    x2 = block2_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,3)
plotDiff(times(xAxis),block1_neutral,block2_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Pleasant Block 1-2')

% Pleasant BLOCK 2-3
for iFrame = 1:size(block1_pleasant,1)
    x1 = block2_neutral(iFrame,:);
    x2 = block3_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,xAxis,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,4)
plotDiff(times(xAxis),block2_neutral,block3_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Pleasant Block 2-3')

%% at subject level

% for iSub = 1:nSub
%     fprintf('Subject %2.2d \n', iSub)
%     tmp = load('-mat', fullfile(STUDY.datasetinfo(iSub).filepath, sprintf('sub-%2.2d.daterp',iSub)));
% 
%     erpData = tmp.chan64;   % O2
%     times = tmp.times;
%     prestim = times > -1000 & times < 10;
%     conds = [tmp.trialinfo.type]';
    
%     figure('color','w')

    % Block 1
%     subplot(3,1,1)

%     for iFrame = 1:size(ERP1,1)
%         x1 = ERP1(iFrame,:);
%         x2 = ERP2(iFrame,:);
%         h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
%         if h1==0 && h2==0 %&& h3==0
% %             cprintf('blue','Normal distribution --> two-sample t-test. \n');
%             [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(2,x1,x2,0.05);
%         else
% %             cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
% %             [p(iFrame),h,stats] = signrank(y,X)
%             [~,~,~,~,~,p(iFrame)] = yuen(x1,x2,10);
%         end
%     end
%     h = p<0.05;
%     plotDiff(times(prestim),ERP1,ERP2,'trimmed mean', h,'Pleasant','Neutral')
%     title('Block 1')
% 
%     % Block 2
%     subplot(3,1,2)
%     idx = blockSize+1:blockSize*2;
%     idxCond1 = strcmp(string(conds(idx,:)), '2');
%     idxCond2 = strcmp(string(conds(idx,:)), '4');
% %     idxCond3 = strcmp(string(conds(idx1,:)), '8');
%     fprintf('Pleasant:    %d trials \n', sum(idxCond1))
%     fprintf('Neutral:     %d trials \n', sum(idxCond2))
% %     fprintf('Unpleasant:  %d trials \n', sum(idxCond3))
%     ERP1 = erpData(prestim,idxCond1);
%     ERP2 = erpData(prestim,idxCond2);
% %     ERP3 = erpData(prestim,idxCond3);
%     for iFrame = 1:size(ERP1,1)
%         x1 = ERP1(iFrame,:);
%         x2 = ERP2(iFrame,:);
%         h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
%         if h1==0 && h2==0 %&& h3==0
% %             cprintf('blue','Normal distribution --> two-sample t-test. \n');
%             [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(2,x1,x2,0.05);
%         else
% %             cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
% %             [p(iFrame),h,stats] = signrank(y,X)
%             [~,~,~,~,~,p(iFrame)] = yuen(x1,x2,10);
%         end
%     end
%     h = p<0.05;
%     plotDiff(times(prestim),ERP1,ERP2,'trimmed mean', h,'Pleasant','Neutral')
%     title('Block 2')
% 
%     % Block 3
%     subplot(3,1,3)
%     idx = blockSize*2+1:blockSize*3;
%     idxCond1 = strcmp(string(conds(idx,:)), '2');
%     idxCond2 = strcmp(string(conds(idx,:)), '4');
% %     idxCond3 = strcmp(string(conds(idx1,:)), '8');
%     fprintf('Pleasant:    %d trials \n', sum(idxCond1))
%     fprintf('Neutral:     %d trials \n', sum(idxCond2))
% %     fprintf('Unpleasant:  %d trials \n', sum(idxCond3))
%     ERP1 = erpData(prestim,idxCond1);
%     ERP2 = erpData(prestim,idxCond2);
% %     ERP3 = erpData(prestim,idxCond3);
%     for iFrame = 1:size(ERP1,1)
%         x1 = ERP1(iFrame,:);
%         x2 = ERP2(iFrame,:);
%         h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
%         if h1==0 && h2==0 %&& h3==0
% %             cprintf('blue','Normal distribution --> two-sample t-test. \n');
%             [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(2,x1,x2,0.05);
%         else
% %             cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
% %             [p(iFrame),h,stats] = signrank(y,X)
%             [~,~,~,~,~,p(iFrame)] = yuen(x1,x2,10);
%         end
%     end
%     h = p<0.05;
%     plotDiff(times(prestim),ERP1,ERP2,'trimmed mean', h,'Pleasant','Neutral')
%     title('Block 3')
%     
% 
% end



%% LIMO files

cd 'D:\presentiment_eeg\data_clean\derivatives\LIMO_presentiment_eeg\pleasant_neutral'
tmp = load('LIMO.mat');
gpLIMO = tmp.LIMO;
eeglab; close;

xAxis = gpLIMO.data.timevect;
peakTime = find(xAxis == -148);

iSub = 1%:78


    load(fullfile(subPath, 'LIMO.mat')) % Structure containing all data and model information
    load(fullfile(subPath, 'Betas.mat')) % Beta parameter estimates
    load(fullfile(subPath, 'Condition_effect_1.mat')) % Factor effect (Categorical designs)
    load(fullfile(subPath, 'Yr.mat')) % Single trial data reorganized to fit X (grouped by condition)
    load(fullfile(subPath, 'Yhat.mat'))  % Predicted data 
    load(fullfile(subPath, 'Res.mat')) % Residuals (non modelled) data 

    conditions = LIMO.data.Cat; % Trial condition
    fvals = Condition_effect(64,:,1);
    pvals = Condition_effect(64,:,2);
    
    sigLat = find(pvals < 0.05);
    xAxis(sigLat)
    
%     figure('color','w'); hold on;
%     for iTrial = 1:size(Yr,3)
%         if conditions(iTrial) == 1
%             plot(t,squeeze(Yr(64,:,iTrial)),'g'); 
%             plot(t,squeeze(Yhat(64,:,iTrial)),'g'); 
%             plot(t,squeeze(Res(64,:,:)),'g'); 
%         elseif conditions(iTrial) == 2
%             plot(t,squeeze(Yr(64,:,iTrial)),'k'); 
%             plot(t,squeeze(Yhat(64,:,iTrial)),'k'); 
%             plot(t,squeeze(Res(64,:,:)),'k'); 
%         elseif conditions(iTrial) == 3
%             plot(t,squeeze(Yr(64,:,iTrial)),'r'); 
%             plot(t,squeeze(Yhat(64,:,iTrial)),'r'); 
%             plot(t,squeeze(Res(64,:,:)),'r'); 
%         end
%     end
    idxCond1 = conditions == 1;
    idxCond2 = conditions == 2;
    idxCond3 = conditions == 3;

    figure('color','w'); 
    subplot(2,2,1)
    plotDiff(xAxis, squeeze(Res(64,:,idxCond1)), squeeze(Res(64,:,idxCond2)), [], 'Pleasant', 'Neutral')
    title('Residuals')
    subplot(2,2,2)
    plotDiff(xAxis, squeeze(Yhat(64,:,idxCond1)), squeeze(Yhat(64,:,idxCond2)), [], 'Pleasant', 'Neutral')
    title('Predicted data')
    subplot(2,2,3)
    plotDiff(xAxis, squeeze(Res(64,:,idxCond3)), squeeze(Res(64,:,idxCond2)), [], 'Unpleasant', 'Neutral')
    title('Residuals')
    subplot(2,2,4)
    plotDiff(xAxis, squeeze(Yhat(64,:,idxCond3)), squeeze(Yhat(64,:,idxCond2)), [], 'Unpleasant', 'Neutral')
    title('Predicted data')
