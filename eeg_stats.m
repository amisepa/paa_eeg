%% Predictive anticipatory activity (PAA) - Stats and results

clear; close all;clc
folder = 'D:\presentiment_eeg\data_clean';
eeglab; close;

% Load study
[STUDY, ALLEEG] = pop_loadstudy('filename','presentiment_eeg.study','filepath',folder);
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% Stats with EEGLAB
STUDY = pop_statparams(STUDY, 'condstats','on','method','perm','mcorrect','fdr','alpha',0.01);
STUDY = pop_erpparams(STUDY, 'plotconditions','together','timerange',[-1950 1950] ,'averagechan','rms');
std_erpplot(STUDY,ALLEEG,'channels',{EEG(1).chanlocs.labels}, 'design', 1);

% Compare conditions with SE
channel = {'Fz'};
[STUDY, erpdata, erptimes] = std_erpplot(STUDY, ALLEEG, 'channels', channel, 'timerange', [-1950 1950], 'noplot', 'on');
std_plotcurve(erptimes, erpdata, 'plotconditions', 'together', 'plotstderr', 'on', 'figure', 'on', 'filter', 10);
legend('Pleasant', 'Neutral', 'Unpleasant'); title([char(channel) ' (mean + SE)']);

% Hierarchichal linear modeling - Compute 1st level (post-stimulus)
pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp', ...
    'timelim',[-50 1950],'erase','on','splitreg','off','interaction','off');

pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp', ...
    'timelim',[-1950 50],'erase','on','splitreg','off','interaction','off');

% 2nd level (post-stimulus)
limo_eeg

