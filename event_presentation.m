%% check events randomness and lags
clear; close all; clc;
codeFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\code';
cleanDataFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\data_clean';
load(fullfile(codeFolder,'sInfo2.mat'));
eeglab
resultFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\events_analysis';
cd(resultFolder);

sInfo(isnan([sInfo.trials])) = [];

for iSub = 1:length(sInfo)
    %raw    
    if sInfo(iSub).group == 1
        dataFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\data_source\group-01';
    else
        dataFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\data_source\group-02';
    end
    EEG = pop_loadcnt(fullfile(dataFolder, [sInfo(iSub).filename '.cnt']), 'dataformat', 'auto', 'memmapfile', '');
    sRate = EEG.srate;
    if ischar(EEG.event(1).type)
        for iEv = 1:length(EEG.event)
            EEG.event(iEv).type = str2double(EEG.event(iEv).type);
        end
    end
    
    count1 = 1;
    count2 = 1;
    count4 = 1;
    count8 = 1;
    for iEv = 2:length(EEG.event)
        if EEG.event(iEv).type == 1 
            cond_prec_stim1(iSub,count1) = EEG.event(iEv-1).type; %condition preceding event 1
            trial1_prevLag(iSub,count1) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between 2 stims in s
            if EEG.event(iEv-1).type == 1
                stim1_lag_stim1(iSub,count1) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 1 and stim 1
            elseif EEG.event(iEv-1).type == 2
                stim2_lag_stim1(iSub,count1) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 2 and stim 1
            elseif EEG.event(iEv-1).type == 4
                stim4_lag_stim1(iSub,count1) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 4 and stim 1
            elseif EEG.event(iEv-1).type == 8
                stim8_lag_stim1(iSub,count1) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 8 and stim 1
            end
            count1 = count1+1;
        elseif EEG.event(iEv).type == 2
            cond_prec_stim2(iSub,count2) = EEG.event(iEv-1).type;
            trial2_prevLag(iSub,count2) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate;
            if EEG.event(iEv-1).type == 1
                stim1_lag_stim2(iSub,count2) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 1 and stim 2
            elseif EEG.event(iEv-1).type == 2
                stim2_lag_stim2(iSub,count2) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 2 and stim 2
            elseif EEG.event(iEv-1).type == 4
                stim4_lag_stim2(iSub,count2) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 4 and stim 2
            elseif EEG.event(iEv-1).type == 8
                stim8_lag_stim2(iSub,count2) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 8 and stim 2
            end
            count2 = count2+1;
        elseif EEG.event(iEv).type == 4 
            cond_prec_stim4(iSub,count4) = EEG.event(iEv-1).type;
            trial4_prevLag(iSub,count4) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate;
            if EEG.event(iEv-1).type == 1
                stim1_lag_stim4(iSub,count4) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 1 and stim 4
            elseif EEG.event(iEv-1).type == 2
                stim2_lag_stim4(iSub,count4) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 2 and stim 4
            elseif EEG.event(iEv-1).type == 4
                stim4_lag_stim4(iSub,count4) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 4 and stim 4
            elseif EEG.event(iEv-1).type == 8
                stim8_lag_stim4(iSub,count4) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 8 and stim 4
            end
            count4 = count4+1;
        elseif EEG.event(iEv).type == 8 
            cond_prec_stim8(iSub,count8) = EEG.event(iEv-1).type;
            trial8_prevLag(iSub,count8) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate;
            if EEG.event(iEv-1).type == 1
                stim1_lag_stim8(iSub,count8) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 1 and stim 8
            elseif EEG.event(iEv-1).type == 2
                stim2_lag_stim8(iSub,count8) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 2 and stim 8
            elseif EEG.event(iEv-1).type == 4
                stim4_lag_stim8(iSub,count8) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 4 and stim 8
            elseif EEG.event(iEv-1).type == 8
                stim8_lag_stim8(iSub,count8) = (EEG.event(iEv).latency - EEG.event(iEv-1).latency) / sRate; %period between stim 8 and stim 8
            end
            count8 = count8+1;
        end
    end
    
%     %event and lags preceding each condition
%     trial2_prev2(iSub,:) = (sum(cond_prec_stim2(iSub,:) == 2) / length(cond_prec_stim2))*100; %trial 2 preceded by event 2 
%     trial2_prev4(iSub,:) = (sum(cond_prec_stim2(iSub,:) == 4) / length(cond_prec_stim2))*100; %trial 2 preceded by event 4 
%     trial2_prev8(iSub,:) = (sum(cond_prec_stim2(iSub,:) == 8) / length(cond_prec_stim2))*100; %trial 2 preceded by event 8 
%     trial2_meanLag(iSub,:) = mean(trial2_prevLag(iSub,:)); 
%     trial2_sdLag(iSub,:) = std(trial2_prevLag(iSub,:));
% 
%     trial4_prev2(iSub,:) = (sum(cond_prec_stim4(iSub,:) == 2) / length(cond_prec_stim4))*100; %trial 4 preceded by event 2 
%     trial4_prev4(iSub,:) = (sum(cond_prec_stim4(iSub,:) == 4) / length(cond_prec_stim4))*100; %trial 4 preceded by event 4 
%     trial4_prev8(iSub,:) = (sum(cond_prec_stim4(iSub,:) == 8) / length(cond_prec_stim4))*100; %trial 4 preceded by event 8 
%     trial4_meanLag(iSub,:) = mean(trial4_prevLag(iSub,:));
%     trial4_sdLag(iSub,:) = std(trial4_prevLag(iSub,:));
% 
%     trial8_prev2(iSub,:) = (sum(cond_prec_stim8(iSub,:) == 2) / length(cond_prec_stim8))*100; %trial 8 preceded by event 2 
%     trial8_prev4(iSub,:) = (sum(cond_prec_stim8(iSub,:) == 4) / length(cond_prec_stim8))*100; %trial 8 preceded by event 4 
%     trial8_prev8(iSub,:) = (sum(cond_prec_stim8(iSub,:) == 8) / length(cond_prec_stim8))*100; %trial 8 preceded by event 8 
%     trial8_meanLag(iSub,:) = mean(trial8_prevLag(iSub,:));
%     trial8_sdLag(iSub,:) = std(trial8_prevLag(iSub,:));
    
    %proportion of clean trials per condition
    subjFolder = fullfile(cleanDataFolder, ['sub-' sInfo(iSub).id], ['ses-' sInfo(iSub).session], 'eeg');
    EEG = pop_loadset('filename', [sInfo(iSub).filename '_eeg.set'],'filepath', subjFolder);
    if ischar(EEG.event(1).type)
        for iEv = 1:length(EEG.event)
            EEG.event(iEv).type = str2double(EEG.event(iEv).type);
        end
    end
    nTrials_all(iSub,:) = sInfo(iSub).trials;
    nTrials_2(iSub,:) = (sum([EEG.event.type] == 2) / nTrials_all(iSub,:)) * 100;
    nTrials_4(iSub,:) = (sum([EEG.event.type] == 4) / nTrials_all(iSub,:)) * 100;
    nTrials_8(iSub,:) = (sum([EEG.event.type] == 8) / nTrials_all(iSub,:)) * 100;
        
end

%% SAVE

%Save number of trials per condition
save(fullfile(resultFolder, 'trial_number', 'nTrials_all.mat'), 'nTrials_all');
writematrix(nTrials_all, fullfile(resultFolder, 'trial_number', 'nTrials_all.xlsx'));
save(fullfile(resultFolder, 'trial_number', 'nTrials_2.mat'), 'nTrials_2');
writematrix(nTrials_2, fullfile(resultFolder, 'trial_number', 'nTrials_2.xlsx'));
save(fullfile(resultFolder, 'trial_number', 'nTrials_4.mat'), 'nTrials_4');
writematrix(nTrials_4, fullfile(resultFolder, 'trial_number', 'nTrials_4.xlsx'));
save(fullfile(resultFolder, 'trial_number', 'nTrials_8.mat'), 'nTrials_8');
writematrix(nTrials_8, fullfile(resultFolder, 'trial_number', 'nTrials_8.xlsx'));

%Save condition preceding each stimulus
save(fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim1.mat'), 'cond_prec_stim1');
writematrix(cond_prec_stim1, fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim1.xlsx'));
save(fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim2.mat'), 'cond_prec_stim2');
writematrix(cond_prec_stim2, fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim2.xlsx'));
save(fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim4.mat'), 'cond_prec_stim4');
writematrix(cond_prec_stim4, fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim4.xlsx'));
save(fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim8.mat'), 'cond_prec_stim8');
writematrix(cond_prec_stim8, fullfile(resultFolder, 'preceding_condition', 'cond_prec_stim8.xlsx'));

%Save period preceding each stimulus
save(fullfile(resultFolder, 'preceding_period', 'period_prec_stim1.mat'), 'trial1_prevLag');
writematrix(trial1_prevLag, fullfile(resultFolder, 'preceding_period', 'period_prec_stim1.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'period_prec_stim2.mat'), 'trial2_prevLag');
writematrix(trial2_prevLag, fullfile(resultFolder, 'preceding_period', 'period_prec_stim2.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'period_prec_stim4.mat'), 'trial4_prevLag');
writematrix(trial4_prevLag, fullfile(resultFolder, 'preceding_period', 'period_prec_stim4.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'period_prec_stim8.mat'), 'trial8_prevLag');
writematrix(trial8_prevLag, fullfile(resultFolder, 'preceding_period', 'period_prec_stim8.xlsx'));

%Save period preceding each stimulus 1 and which stim it was
save(fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim1.mat'), 'stim1_lag_stim1');
writematrix(stim1_lag_stim1, fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim1.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim1.mat'), 'stim2_lag_stim1');
writematrix(stim2_lag_stim1, fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim1.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim1.mat'), 'stim4_lag_stim1');
writematrix(stim4_lag_stim1, fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim1.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim1.mat'), 'stim8_lag_stim1');
writematrix(stim8_lag_stim1, fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim1.xlsx'));

%Save period preceding each stimulus 2 and which stim it was
save(fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim2.mat'), 'stim1_lag_stim2');
writematrix(stim1_lag_stim2, fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim2.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim2.mat'), 'stim2_lag_stim2');
writematrix(stim2_lag_stim2, fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim2.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim2.mat'), 'stim4_lag_stim2');
writematrix(stim4_lag_stim2, fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim2.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim2.mat'), 'stim8_lag_stim2');
writematrix(stim8_lag_stim2, fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim2.xlsx'));

%Save period preceding each stimulus 4 and which stim it was
save(fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim4.mat'), 'stim1_lag_stim4');
writematrix(stim1_lag_stim4, fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim4.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim4.mat'), 'stim2_lag_stim4');
writematrix(stim2_lag_stim4, fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim4.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim4.mat'), 'stim4_lag_stim4');
writematrix(stim4_lag_stim4, fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim4.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim4.mat'), 'stim8_lag_stim4');
writematrix(stim8_lag_stim4, fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim4.xlsx'));

%Save period preceding each stimulus 8 and which stim it was
save(fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim8.mat'), 'stim1_lag_stim8');
writematrix(stim1_lag_stim8, fullfile(resultFolder, 'preceding_period', 'stim1_lag_stim8.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim8.mat'), 'stim2_lag_stim8');
writematrix(stim2_lag_stim8, fullfile(resultFolder, 'preceding_period', 'stim2_lag_stim8.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim8.mat'), 'stim4_lag_stim8');
writematrix(stim4_lag_stim8, fullfile(resultFolder, 'preceding_period', 'stim4_lag_stim8.xlsx'));
save(fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim8.mat'), 'stim8_lag_stim8');
writematrix(stim8_lag_stim8, fullfile(resultFolder, 'preceding_period', 'stim8_lag_stim8.xlsx'));

%% Number of clean trials (total and per condition)

%ALL conditions
nTrials_all_mean = round(mean(nTrials_all));
nTrials_all_sd = round(std(nTrials_all));
figure; set(gcf,'Color','w');
h = histogram(nTrials_all); % 'BinWidth',1 %'Orientation', 'horizontal'
title(['Number of clean trials per subject (mean ' num2str(nTrials_all_mean) ' +- ' num2str(nTrials_all_sd) ')']);
ylabel('Number of trials'); xlabel('Number of subjects');
E = h.BinEdges; y = h.BinCounts; xloc = E(1:end-1)+diff(E)/2; text(xloc, y+1, string(y)) %display values above bars
h.FaceColor = [0 0.5 0.5]; 

%Pleasant
nTrials_2_mean = round(mean(nTrials_2));
nTrials_2_sd = round(std(nTrials_2));
figure; set(gcf,'Color','w');
h = histogram(nTrials_2); % 'BinWidth',1
title(['Number of clean Pleasant trials per subject (mean ' num2str(nTrials_2_mean) ' +- ' num2str(nTrials_2_sd) ')']);
h.FaceColor = [0 0.5 0.5]; 
E = h.BinEdges;y = h.BinCounts;xloc = E(1:end-1)+diff(E)/2;text(xloc, y+1, string(y)) %display values above bars
ylabel('Number of clean trials'); xlabel('Number of subjects');

%Neutral
nTrials_4_mean = round(mean(nTrials_4));
nTrials_4_sd = round(std(nTrials_4));
figure; set(gcf,'Color','w');
h = histogram(nTrials_4); % 'BinWidth',1
title(['Number of clean Neutral trials per subject (mean ' num2str(nTrials_4_mean) ' +- ' num2str(nTrials_4_sd) ')']);
h.FaceColor = [0 0.5 0.5];
E = h.BinEdges;y = h.BinCounts;xloc = E(1:end-1)+diff(E)/2;text(xloc, y+1, string(y)) %display values above bars
ylabel('Number of clean trials'); xlabel('Number of subjects');

%Unpleasant
nTrials_8_mean = round(mean(nTrials_8));
nTrials_8_sd = round(std(nTrials_8));
figure; set(gcf,'Color','w');
h = histogram(nTrials_8); % 'BinWidth',1
title(['Number of clean Unpleasant trials per subject (mean ' num2str(nTrials_8_mean) ' +- ' num2str(nTrials_8_sd) ')']);
h.FaceColor = [0 0.5 0.5];
E = h.BinEdges;y = h.BinCounts;xloc = E(1:end-1)+diff(E)/2;text(xloc, y+1, string(y)) %display values above bars
ylabel('Number of clean trials'); xlabel('Number of subjects');

%% Condition preceding each stimulus
% histogram probability option: The height of each bar is the relative number 
% of observations (number of observations  in bin / total number of
% observations). The sum of the bar heights is less than or equal to 1.

cond_prec_stim2((cond_prec_stim2 == 0 | cond_prec_stim2 == 255 | cond_prec_stim2 == 223 | cond_prec_stim2 == 239 | cond_prec_stim2 == 207)) = NaN;
prev_2 = categorical(cond_prec_stim2,[1 2 4 8], {'Checkerboard','Pleasant','Neutral','Unpleasant'});
figure; set(gcf,'Color','w');
h = histogram(prev_2, 'BarWidth', 0.5, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; ylabel('%');
title('Events preceding Pleasant stimuli (all trials and all subjects)');

cond_prec_stim4((cond_prec_stim4 == 0 | cond_prec_stim4 == 255 | cond_prec_stim4 == 223 | cond_prec_stim4 == 239 | cond_prec_stim4 == 207)) = NaN;
prev_4 = categorical(cond_prec_stim4,[1 2 4 8], {'Checkerboard','Pleasant','Neutral','Unpleasant'});
figure; set(gcf,'Color','w');
h = histogram(prev_4, 'BarWidth', 0.5, 'Normalization', 'probability'); 
title('Events preceding Neutral stimuli (all trials and all subjects)');
h.FaceColor = [0 0.5 0.5]; xlabel('%');

cond_prec_stim8((cond_prec_stim8 == 0 | cond_prec_stim8 == 255 | cond_prec_stim8 == 223 | cond_prec_stim8 == 239 | cond_prec_stim8 == 207)) = NaN;
prev_8 = categorical(cond_prec_stim8,[1 2 4 8], {'Checkerboard','Pleasant','Neutral','Unpleasant'});
figure; set(gcf,'Color','w');
h = histogram(prev_8, 'BarWidth', 0.5, 'Normalization', 'probability');
title('Events preceding Unpleasant stimuli (all trials and all subjects)');
h.FaceColor = [0 0.5 0.5]; xlabel('%');

%% Period preceding each stimulus (in s)

data = trial1_prevLag;
name = 'period_prec_stim1';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period preceding Checkerboard stimulus (all trials and all subjects)');
% E = h.BinEdges; y = h.Values; xloc = E(1:end-1)+diff(E)/2; text(xloc, y+1, string(y)) %display values above bars
xlim([min(min(data)) max(max(data))]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = trial2_prevLag;
name = 'period_prec_stim2';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period preceding Pleasant stimulus (all trials and all subjects)');
% E = h.BinEdges; y = h.Values; xloc = E(1:end-1)+diff(E)/2; text(xloc, y+1, string(y)) %display values above bars
xlim([min(min(data)) max(max(data))]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = trial4_prevLag;
name = 'period_prec_stim4';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period preceding Neutral stimulus (all trials and all subjects)');
% E = h.BinEdges; y = h.Values; xloc = E(1:end-1)+diff(E)/2; text(xloc, y+1, string(y)) %display values above bars
xlim([min(min(data)) max(max(data))]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = trial8_prevLag;
name = 'period_prec_stim8';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period preceding Unpleasant stimulus (all trials and all subjects)');
% E = h.BinEdges; y = h.Values; xloc = E(1:end-1)+diff(E)/2; text(xloc, y+1, string(y)) %display values above bars
xlim([min(min(data)) max(max(data))]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

%% Period preceding each stimulus for each condition (in s)

data = stim1_lag_stim2;
name = 'stim1_lag_stim2';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Checkerboard and Pleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim2_lag_stim2;
name = 'stim2_lag_stim2';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Pleasant and Pleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim4_lag_stim2;
name = 'stim4_lag_stim2';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Neutral and Pleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim8_lag_stim2;
name = 'stim8_lag_stim2';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Unplesant and Pleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim1_lag_stim4;
name = 'stim1_lag_stim4';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Checkerboard and Neutral stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim2_lag_stim4;
name = 'stim2_lag_stim4';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Pleasant and Neutral stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim4_lag_stim4;
name = 'stim4_lag_stim4';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Neutral and Neutral stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim8_lag_stim4;
name = 'stim8_lag_stim4';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Unpleasant and Neutral stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim1_lag_stim8;
name = 'stim1_lag_stim8';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Checkerboard and Unpleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim2_lag_stim8;
name = 'stim2_lag_stim8';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Pleasant and Unpleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim4_lag_stim8;
name = 'stim4_lag_stim8';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Neutral and Unpleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

data = stim8_lag_stim8;
name = 'stim8_lag_stim8';
figure; set(gcf,'Color','w');
data(data == 0) = NaN;
h = histogram(data(~isnan(data)), 'BinWidth',0.1, 'Normalization', 'probability');
h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
title('Period between Unpleasant and Unpleasant stim');
xlim([min(min(data)) max(max(data))]);ylim([0 1]);
saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));

% data = stim1_lag_stim1;
% name = 'stim1_lag_stim1';
% figure; set(gcf,'Color','w');
% h = histogram(data, 'BinWidth',1, 'Normalization', 'probability');
% h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
% title('Period between Checkerboard stim and previous checkerboard stim');
% xlim([min(min(data)) max(max(data))]);ylim([0 1]);
% saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
% print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));
% 
% data = stim2_lag_stim1;
% name = 'stim2_lag_stim1';
% figure; set(gcf,'Color','w');
% h = histogram(data, 'BinWidth',1, 'Normalization', 'probability');
% h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
% title('Period between Pleasant stim and previous Checkerboard stim');
% xlim([min(min(data)) max(max(data))]); ylim([0 1]);
% saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
% print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));
% 
% data = stim4_lag_stim1;
% name = 'stim4_lag_stim1';
% figure; set(gcf,'Color','w');
% h = histogram(data, 'BinWidth',1, 'Normalization', 'probability');
% h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
% title('Period between Neutral stim and previous Checkerboard stim');
% xlim([min(min(data)) max(max(data))]); ylim([0 1]);
% saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
% print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));
% 
% data = stim8_lag_stim1;
% name = 'stim8_lag_stim1';
% figure; set(gcf,'Color','w');
% h = histogram(data, 'BinWidth',1, 'Normalization', 'probability');
% h.FaceColor = [0 0.5 0.5]; xlabel('Time (s)'); ylabel('%');
% title('Period between Unpleasant stim and previous Checkerboard stim');
% xlim([min(min(data)) max(max(data))]); ylim([0 1]);
% saveas(h, fullfile(resultFolder, 'preceding_period', [name '.fig']));
% print(gcf, '-dtiffn', fullfile(resultFolder, 'preceding_period', [name '.tif']));


