%% Check time-on-task effects (i.e., how effect may change over experiment)

clear; close all;clc
dataDir = 'D:\presentiment_eeg\data_clean';
outDir = 'D:\presentiment_eeg\order_effects';
cd(dataDir)
eeglab; close;

% Load study
[STUDY, ALLEEG] = pop_loadstudy('filename','presentiment_eeg.study','filepath',dataDir);

% Calculate minmimum block size
tmp = { STUDY.datasetinfo.trialinfo }';
nSub = length(tmp);
for i = 1:nSub
    trialLength(i,:) = length(STUDY.datasetinfo(i).trialinfo);
end
minTrial = min(trialLength);
blockSize = floor(minTrial / 3);
clear tmp

%% Get ERP data for 3 blocks of trials

% chan = 64;   % 64 = O2

% [STUDY, erpdata, erptimes] = std_erpplot(STUDY,ALLEEG,'channels',{EEG.chanlocs(chan).labels},'timerange',[-1000 10],'noplot','on');

for iSub = 1:nSub
    fprintf('Subject %2.2d \n', iSub)
    tmp = load('-mat', fullfile(STUDY.datasetinfo(iSub).filepath, sprintf('sub-%2.2d.daterp',iSub)));
        
    times = tmp.times;
    prestim = times > -1000 & times < 1000;
    
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

    % block 2
    idx = blockSize+1:blockSize*2;
    idxCond1 = strcmp(string(conds(idx,:)), '2');
    idxCond2 = strcmp(string(conds(idx,:)), '4');
    idxCond3 = strcmp(string(conds(idx,:)), '8');
    block2_nTrials(iSub,1) = sum(idxCond1);
    block2_nTrials(iSub,2) = sum(idxCond2);
    block2_nTrials(iSub,3) = sum(idxCond3);

    % block 3
    idx = blockSize*2+1:blockSize*3;
    idxCond1 = strcmp(string(conds(idx,:)), '2');
    idxCond2 = strcmp(string(conds(idx,:)), '4');
    idxCond3 = strcmp(string(conds(idx,:)), '8');
    block3_nTrials(iSub,1) = sum(idxCond1);
    block3_nTrials(iSub,2) = sum(idxCond2);
    block3_nTrials(iSub,3) = sum(idxCond3);
    
    % ERP data for each condition and block
    for iChan = 1:64
    
        erpData = tmp{:,iChan};   
            
        % block 1
        block1_pleasant(iChan,:,iSub) = trimmean(erpData(prestim,idxCond1),10,2);
        block1_neutral(iChan,:,iSub) = trimmean(erpData(prestim,idxCond2),10,2);
        block1_unplesant(iChan,:,iSub) = trimmean(erpData(prestim,idxCond3),10,2);

        % block 2
        block2_pleasant(iChan,:,iSub) = trimmean(erpData(prestim,idxCond1),10,2);
        block2_neutral(iChan,:,iSub) = trimmean(erpData(prestim,idxCond2),10,2);
        block2_unplesant(iChan,:,iSub) = trimmean(erpData(prestim,idxCond3),10,2);
    
        % block 3
        block3_pleasant(iChan,:,iSub) = trimmean(erpData(prestim,idxCond1),10,2);
        block3_neutral(iChan,:,iSub) = trimmean(erpData(prestim,idxCond2),10,2);
        block3_unplesant(iChan,:,iSub) = trimmean(erpData(prestim,idxCond3),10,2);
    end

end

% % Compare number of trials per condition per block
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

disp('done')

%% Compare conditions for each block

mcctype = 2;
chanlocs = ALLEEG(1).chanlocs;

% Block 1
[results, results_H0] = compute_randomeffect(block1_pleasant,block1_neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'block1_pairwise_comparisons.mat'), 'results')
save(fullfile(outDir, 'block1_pairwise_comparisons_H0.mat'), 'results_H0', '-v7.3')
[neighbors, neighbormatrix] = get_channelneighbors(chanlocs,0);
t1 = tic;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
fprintf('Computing time: %f min \n', toc(t1)/60)
if sum(mask) > 0
    plotresults('time', times(prestim), tvals, mask, pcorr, mcctype, chanlocs)
end

% Block 2
[results, results_H0] = compute_randomeffect(block2_pleasant,block2_neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'block2_pairwise_comparisons.mat'), 'results')
save(fullfile(mainfolder,'block2_pairwise_comparisons_H0.mat'), 'results_H0', '-v7.3')
t1 = tic;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
fprintf('Computing time: %f min \n', toc(t1)/60)
if sum(mask) > 0
    plotresults('time', times(prestim), tvals, mask, pcorr, mcctype, chanlocs)
end

% Block 3
[results, results_H0] = compute_randomeffect(block3_pleasant,block3_neutral,1000,'trimmed mean');
tvals = results(:,:,1);
pvals = results(:,:,2);
tvals_H0 = squeeze(results_H0(:,:,1,:));
pvals_H0 = squeeze(results_H0(:,:,2,:));
save(fullfile(outDir, 'block3_pairwise_comparisons.mat'), 'results')
save(fullfile(mainfolder,'block3_pairwise_comparisons_H0.mat'), 'results_H0', '-v7.3')
t1 = tic;
[mask, pcorr, nClust] = compute_mcc(tvals, pvals, tvals_H0, pvals_H0, mcctype, 0.05, neighbormatrix);
fprintf('Computing time: %f min \n', toc(t1)/60)
if sum(mask) > 0
    plotresults('time', times(prestim), tvals, mask, pcorr, mcctype, chanlocs)
end

gong

%% Stats comparing across conditions for each block

% ERP BLOCK 1
figure('color','w');
for iFrame = 1:size(block1_pleasant,1)
    x1 = block1_pleasant(iFrame,:);
    x2 = block1_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        % cprintf('blue','Normal distribution --> two-sample t-test. \n');
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(3,1,1)
plotDiff(times(prestim),block1_pleasant,block1_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Block 1')

% ERP BLOCK 2
for iFrame = 1:size(block2_pleasant,1)
    x1 = block2_pleasant(iFrame,:);
    x2 = block2_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        % cprintf('blue','Normal distribution --> two-sample t-test. \n');
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(3,1,2)
plotDiff(times(prestim),block2_pleasant,block2_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Block 2')

% ERP BLOCK 3
for iFrame = 1:size(block3_pleasant,1)
    x1 = block3_pleasant(iFrame,:);
    x2 = block3_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        % cprintf('blue','Normal distribution --> two-sample t-test. \n');
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  cprintf('blue','Non-normal distribution --> Yuen t-test. \n');
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(3,1,3)
plotDiff(times(prestim),block3_pleasant,block3_neutral,'trimmed mean', h,'Pleasant','Neutral')
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
    prestim = times > -1000 & times < 1000;
    conds = [tmp.trialinfo.type]';

    idxCond1 = strcmp(string(conds(idx,:)), '2');
    idxCond2 = strcmp(string(conds(idx,:)), '4');
    idxCond3 = strcmp(string(conds(idx,:)), '8');

    pleasant(:,iSub) = trimmean(erpData(prestim,idxCond1),20,2);
    neutral(:,iSub) = trimmean(erpData(prestim,idxCond2),20,2);
    unpleasant(:,iSub) = trimmean(erpData(prestim,idxCond3),20,2);

end

% Pleasant-neutral
figure('color','w');
subplot(2,1,1)
for iFrame = 1:size(pleasant,1)
    x1 = pleasant(iFrame,:);
    x2 = neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2);
    if h1==0 && h2==0
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
h = fdr_bh(p,0.05,'pdep','no');
plotDiff(times(prestim),pleasant,neutral,'trimmed mean', h,'Pleasant','Neutral')
title('All trials (channel O2)')

% Pleasant-neutral
subplot(2,1,2)
for iFrame = 1:size(pleasant,1)
    x1 = unpleasant(iFrame,:);
    x2 = neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2);
    if h1==0 && h2==0
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
h = fdr_bh(p,0.05,'pdep','no');
plotDiff(times(prestim),unpleasant,neutral,'trimmed mean', h,'Unpleasant','Neutral')
title('All trials (channel O2)')

%% Stats comparing same condition across blocks

% Pleasant BLOCK 1-2
figure('color','w');
for iFrame = 1:size(block1_pleasant,1)
    x1 = block1_pleasant(iFrame,:);
    x2 = block2_pleasant(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,1)
plotDiff(times(prestim),block1_pleasant,block2_pleasant,'trimmed mean', h,'Pleasant','Neutral')
title('Pleasant Block 1-2')

% Pleasant BLOCK 2-3
for iFrame = 1:size(block1_pleasant,1)
    x1 = block2_pleasant(iFrame,:);
    x2 = block3_pleasant(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,2)
plotDiff(times(prestim),block2_pleasant,block3_pleasant,'trimmed mean', h,'Pleasant','Neutral')
title('Pleasant Block 2-3')

% Neutral BLOCK 1-2
for iFrame = 1:size(block1_pleasant,1)
    x1 = block1_neutral(iFrame,:);
    x2 = block2_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,3)
plotDiff(times(prestim),block1_neutral,block2_neutral,'trimmed mean', h,'Pleasant','Neutral')
title('Pleasant Block 1-2')

% Pleasant BLOCK 2-3
for iFrame = 1:size(block1_pleasant,1)
    x1 = block2_neutral(iFrame,:);
    x2 = block3_neutral(iFrame,:);
    h1 = lillietest(x1); h2 = lillietest(x2); %h3 = variance_homogeneity(x1,x2);
    if h1==0 && h2==0 %&& h3==0
        [m,dfe,ci,sd,n,t,p(iFrame)] = limo_ttest(1,x1,x2,0.05);
    else
        %  [p(iFrame),h,stats] = signrank(y,X)
        [~,~,~,~,~,p(iFrame)] = yuend(x1,x2,20);
    end
end
% h = p<0.05;
h = fdr_bh(p,0.05,'pdep','no');
subplot(2,2,4)
plotDiff(times(prestim),block2_neutral,block3_neutral,'trimmed mean', h,'Pleasant','Neutral')
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

t = gpLIMO.data.timevect;
peakTime = find(t == -148);

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
    t(sigLat)
    
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
    plotDiff(t, squeeze(Res(64,:,idxCond1)), squeeze(Res(64,:,idxCond2)), [], 'Pleasant', 'Neutral')
    title('Residuals')
    subplot(2,2,2)
    plotDiff(t, squeeze(Yhat(64,:,idxCond1)), squeeze(Yhat(64,:,idxCond2)), [], 'Pleasant', 'Neutral')
    title('Predicted data')
    subplot(2,2,3)
    plotDiff(t, squeeze(Res(64,:,idxCond3)), squeeze(Res(64,:,idxCond2)), [], 'Unpleasant', 'Neutral')
    title('Residuals')
    subplot(2,2,4)
    plotDiff(t, squeeze(Yhat(64,:,idxCond3)), squeeze(Yhat(64,:,idxCond2)), [], 'Unpleasant', 'Neutral')
    title('Predicted data')
