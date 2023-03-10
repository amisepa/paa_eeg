%% Predictive anticipatory activity (PAA) - Stats and results

clear; close all;clc
folder = 'D:\presentiment_eeg\data_clean';
cd(folder)
eeglab; close;

% Load study
[STUDY, ALLEEG] = pop_loadstudy('filename','presentiment_eeg.study','filepath',folder);
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% Stats with EEGLAB
STUDY = pop_statparams(STUDY, 'condstats','on','method','perm','mcorrect','fdr','alpha',0.01);
STUDY = pop_erpparams(STUDY, 'plotconditions','together','timerange',[-1300 1950] ,'averagechan','rms');
% std_erpplot(STUDY,ALLEEG,'channels',{EEG(1).chanlocs.labels}, 'design', 1);

% Compare conditions with SE
channel = {'O2'}; 
[STUDY, erpdata, erptimes] = std_erpplot(STUDY,ALLEEG,'channels',channel,'timerange',[-800 10],'noplot','on');
std_plotcurve(erptimes, erpdata, 'plotconditions', 'together', 'plotstderr', 'on', 'figure', 'on', 'filter', 10);
legend('Pleasant', 'Neutral', 'Unpleasant'); title([char(channel) ' (mean + SE)']);

channel = {'FCz'};
[STUDY, erpdata, erptimes] = std_erpplot(STUDY, ALLEEG, 'channels', channel, 'timerange', [-1300 1950], 'noplot', 'on');
std_plotcurve(erptimes, erpdata, 'plotconditions', 'together', 'plotstderr', 'on', 'figure', 'on', 'filter', 10);
legend('Pleasant', 'Neutral', 'Unpleasant'); title([char(channel) ' (mean + SE)']);


%% Hierarchichal linear modeling - Compute 1st level (pre-stimulus)
pop_limo(STUDY, ALLEEG,'method','WLS','measure','daterp', ...
    'timelim',[-1300 15],'erase','on','splitreg','off','interaction','off');

%% Hierarchichal linear modeling - Compute 1st level (post-stimulus)
% pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp', ...
%     'timelim',[-50 1950],'erase','on','splitreg','off','interaction','off');

%% 2nd level (post-stimulus)

limo_eeg

%% Adjust figure for paper

set(gcf,'Color','White')
set(gca,'FontSize',12,'FontWeight','Bold');

%% Main

title('Pleasant vs. Neutral (uncorrected)')

%% Topo

title('Topography at -148 ms')

%% Course plot

title('T-values at channel O2')

%% Save

figName = inputdlg('Figure name for saving:');

print(gcf, [figName{:} '.png'],'-dpng','-r300');   % 300 dpi .png
print(gcf, [figName{:} '.tiff'],'-dtiff','-r300');  % 300 dpi .tiff

