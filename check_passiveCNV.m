%% Check for passive CNV (time-on-task) effects. 
% Load each file and truncate in 3 blocks. Each block has a session number 
% to run LIMO and compare each block.

clear; close all; clc;
mainDir = 'D:\presentiment_eeg';
dataDir = fullfile(mainDir, 'data_clean');
outDir = 'D:\presentiment_eeg\data_truncated';
cd(outDir)

eeglab; close;

for iSub = 1:78

    disp('--------------------------------------------')
    fprintf('                SUBJECT %g \n', iSub)
    disp('--------------------------------------------')

    % Load file
    filepath = fullfile(dataDir, sprintf('sub-%2.2d',iSub)); 
    filename = [sprintf('sub-%2.2d_eeg',iSub) '.set'];
    tmpeeg = pop_loadset('filename',filename,'filepath',filepath);
    sInfo(iSub).nTrials = size(tmpeeg.data,3);

    % Divide file in 3 blocks
    nTrials = size(tmpeeg.data,3);
    blockSize = floor(nTrials/3);
    EEG1 = pop_rejepoch(tmpeeg, blockSize+1:nTrials, 0);
    EEG2 = pop_rejepoch(tmpeeg, 1:blockSize, 0);
    EEG2 = pop_rejepoch(EEG2, blockSize+1:size(EEG2.data,3), 0);
    EEG3 = pop_rejepoch(tmpeeg, 1:blockSize*2+2, 0);
    
    % Check they have the same size
    fprintf('Number of trials per block after epoching: %g %g %g \n', size(EEG1.data,3), size(EEG2.data,3), size(EEG3.data,3))

%     % Plot ERP
%     figure('color','w'); 
%     subplot(3,1,1)
%     pop_timtopo(EEG1, [-1300 1950], NaN);
%     subplot(3,1,2)
%     pop_timtopo(EEG2, [-1300 1950], NaN);
%     subplot(3,1,3)
%     pop_timtopo(EEG3, [-1300 1950], NaN);
    
    % Save with session names
    newpath = fullfile(outDir, sprintf('sub-%2.2d',iSub), 'ses-01'); mkdir(newpath);
    newname = sprintf('sub-%2.2d_ses-01.set',iSub);
    pop_saveset(EEG1,'filepath',newpath,'filename',newname);

    newpath = fullfile(outDir, sprintf('sub-%2.2d',iSub), 'ses-02'); mkdir(newpath);
    newname = sprintf('sub-%2.2d_ses-02.set',iSub);
    pop_saveset(EEG2,'filepath',newpath,'filename',newname);

    newpath = fullfile(outDir, sprintf('sub-%2.2d',iSub), 'ses-03'); mkdir(newpath);
    newname = sprintf('sub-%2.2d_ses-03.set',iSub);
    pop_saveset(EEG3,'filepath',newpath,'filename',newname);
    
end
gong

%% Create study
clear; close all; clc;
dataDir = 'D:\presentiment_eeg\data_truncated';
eeglab; close;

count = 0;
commands = {};
for iSub = 1:78

    for iSes = 1:3
        filepath = fullfile(dataDir, sprintf('sub-%2.2d',iSub), sprintf('ses-%2.2d',iSes)); 
        filename = sprintf('sub-%2.2d_ses-%2.2d.set',iSub,iSes);
        EEG = pop_loadset('filename',filename,'filepath',filepath);
        EEG.subject = sprintf('sub-%2.2d',iSub); 
        EEG.saved = 'no';
        count = count + 1;
        commands = [ commands(:)' 'index' count 'load' fullfile(filepath, filename) 'session' iSes];
        [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','presentiment_eeg_truncated', ...
            'commands',commands,'updatedat','on','savedat','off','rmclust','off');
        [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); 
        CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
    end
end

[STUDY, ALLEEG] = pop_savestudy(STUDY,EEG,'filename','presentiment_eeg_truncated.study','filepath',dataDir);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
eeglab redraw
gong

%% Study design

STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','STUDY.design 1','delfiles','off', ...
    'defaultdesign','off','variable1','type','values1',{'2','4','8'}, ...
    'vartype1','categorical','variable2','session','values2',{1,2,3}, ...
    'vartype2','continuous','subjselect',{'sub-01'});

% Precompute ERP
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','recompute','on','erp','on');
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% Save
[STUDY,EEG] = pop_savestudy(STUDY,EEG,'filename','presentiment_eeg.study','filepath',dataDir);

%% Hierarchichal linear modeling - Compute 1st level with SPLIT REGRESSION ON

pop_limo(STUDY,ALLEEG,'method','WLS','measure','daterp', ...
    'timelim',[-1300 15],'erase','on','splitreg','on','interaction','off');

%% Hierarchichal linear modeling - Compute 1st level (post-stimulus)
% pop_limo(STUDY, ALLEEG, 'method','WLS','measure','daterp', ...
%     'timelim',[-50 1950],'erase','on','splitreg','off','interaction','off');

%% 2nd level (MANUAL)

limo_eeg

%% Save corrected results

set(gcf,'Color','White')
set(gca,'FontSize',12,'FontWeight','Bold');
figName = inputdlg('Figure name for saving:');
print(gcf, [figName{:} '.png'],'-dpng','-r300');    % .png
print(gcf, [figName{:} '.tiff'],'-dtiff','-r300');  % .tiff 

%% %%%%%%%%%%%%%%%%%%%%%%% CHECK INTER-STIMULUS INTERVAL %%%%%%%%%%%%%%%%%
% We observe an passive CNV effect in block 3 Pleasant-Neutral but not in 
% Unpleasant-Neutral, suggesting our effect is due to expectation (brain
% learning temporal structure of sequence and by chance, pleasant trials
% may appear at highest likelihood in the ISI). 
% To check for this, extract pleasant trials preceded by several neutral
% trials and create a new session file. Do the same with remaining pleasant
% trials, run LIMO on each and compare. 

clear; close all; clc;
mainDir = 'D:\presentiment_eeg';
dataDir = fullfile(mainDir, 'data_raw');
outDir = 'D:\presentiment_eeg\data_truncated';
cd(outDir)
codeDir = fullfile(mainDir, 'paa_eeg');
load(fullfile(codeDir, 'sInfo_old2.mat'));
eeglab; close;

% number of preceding trials to examine
ntrials = 2;  

pleasant = [];
neutral = [];
unpleasant = [];
% emotional = [];
for iSub = 1:78
    
    disp('--------------------------------------------')
    fprintf('                SUBJECT %g \n', iSub)
    disp('--------------------------------------------')

    % Load file
%     filepath = fullfile(dataDir, sprintf('sub-%2.2d',iSub)); 
%     filename = [sprintf('sub-%2.2d_eeg',iSub) '.set'];
%     tmpeeg = pop_loadset('filename',filename,'filepath',filepath);
    if sInfo(iSub).group == 1
        filepath = fullfile(dataDir, 'group-01');
    else
        filepath = fullfile(dataDir, 'group-02');
    end
    EEG = pop_loadcnt(fullfile(filepath, [sInfo(iSub).filename '.cnt']), ...
        'dataformat','auto','keystroke','on','memmapfile','');
    EEG = pop_epoch(EEG, {'2' '4' '8'}, [-1.5 2],'epochinfo', 'yes'); %2 = pleasant pictures; 4 = neutral; 8 = unpleasant; 1 = checkerboard


    % Extract block 3 as before
    nTrials = size(EEG.data,3);
    sInfo(iSub).nTrials = nTrials;
    blockSize = floor(nTrials/3);
    EEG = pop_rejepoch(EEG, 1:blockSize*2-1, 0); % - 1 here includes 3 trials before beginning of block 3 start analysis at trial 1 of block 3

    % Remove weird events
    events = str2double({EEG.event.type});
%     events(isnan(events)) = [];
%     events(events == 255 | events == 239 | events == 223 | events == 207) = [];
%     summary(categorical(events))

    % check randomness of sequence
    [~,p_random(iSub)] = runstest(events);
    if p_random(iSub) < 0.05, warning('This trial sequence is no random'); end
    
    for iEv = 4:length(events)
        % Pleasant condition
        if events(iEv) == 2
            % save preceding trials for each trial            
            Pleasant(iSub,iEv,:) = events(iEv-ntrials:iEv-1); 
            
            % save all preceding trials in one array
            pleasant = [ pleasant events(iEv-ntrials:iEv-1) ]; 

            % track Pleasant trials preceded by 3 neutral ones
            if sum( events(iEv-ntrials:iEv-1) == 4 ) == ntrials
                expectPleasant(iSub,iEv) = 1;
                EEG.event(iEv).type = 'EP';     % expected pleasant 
            else
                expectPleasant(iSub,iEv) = 0;
                EEG.event(iEv).type = 'UP';     % unexpected pleasant
            end

        % Neutral condition
        elseif events(iEv) == 4
            % save preceding trials for each trial
            Neutral(iSub,iEv,:) = events(iEv-ntrials:iEv-1); 

            % save all preceding trials in one array
            neutral = [ neutral events(iEv-ntrials:iEv-1) ];

            % track neutral trials with 3 previous emotional trials
            if sum( events(iEv-ntrials:iEv-1) == 2 | events(iEv-ntrials:iEv-1) == 8 ) == ntrials
                expectNeutral(iSub,iEv) = 1;
                EEG.event(iEv).type = 'EN';     % expected neutral 
            else
                expectNeutral(iSub,iEv) = 0;
                EEG.event(iEv).type = 'UN';     % unexpected neutral 
            end
        
        % Unpleasant condition
        elseif events(iEv) == 8
            % save preceding trials for each trial
            Unpleasant(iSub,iEv,:) = events(iEv-ntrials:iEv-1); 

            % save all preceding trials in one array
            unpleasant = [ unpleasant events(iEv-ntrials:iEv-1) ];

            % track Unpleasant trials preceded by 3 neutral ones
            if sum( events(iEv-ntrials:iEv-1) == 4 ) == ntrials
                expectUnpleasant(iSub,iEv) = 1;
                EEG.event(iEv).type = 'EU';     % expected unpleasant 
            else
                expectUnpleasant(iSub,iEv) = 0;
                EEG.event(iEv).type = 'UU';     % unexpected unpleasant 
            end

%         % Emotional condition
%         elseif events(iEv) == 2 || events(iEv) == 8
% 
%             % track emotional trials preceded by 3 neutral trials
%             emotional =  [ emotional events(iEv-ntrials:iEv-1) ];
%             if sum( events(iEv-ntrials:iEv-1) == 4 ) == ntrials
%                 expectEmotional(iSub,iEv) = 1;
%             else
%                 expectEmotional(iSub,iEv) = 0;
%             end
        end
    end
    
    % remove first 3 trials from block 2 and seperate into 2 files
    % (expected and unexpected)
    EEG = pop_rejepoch(EEG, 1:3, 0); % - 1 here includes 3 trials before beginning of block 3 start analysis at trial 1 of block 3
    pop_eegplot(EEG,1,1,1);
    
    % expected 
    idx = strcmp({EEG.event.type}, 'UP');
    EEG = pop_rejepoch(EEG, 1:3, 0); % - 1 here includes 3 trials before beginning of block 3 start analysis at trial 1 of block 3


    fprintf('%g pleasant trials preceded by %g neutral ones \n', sum(expectPleasant(iSub,:)), ntrials)
    fprintf('%g neutral trials preceded by %g emotional ones \n', sum(expectNeutral(iSub,:)), ntrials)
    fprintf('%g pleasant trials preceded by %g neutral ones \n', sum(expectUnpleasant(iSub,:)), ntrials)
end

figure('color','w');
boxplot([ sum(expectPleasant,2)  sum(expectNeutral,2) sum(expectUnpleasant,2)])
xticklabels({'Pleasant' 'Neutral' 'Unpleasant'})
title('preceding trials for each condition per subject')

% % Remove other events
% pleasant(isnan(pleasant)) = []; neutral(isnan(neutral)) = []; unpleasant(isnan(unpleasant)) = [];
% pleasant(pleasant == 255 | pleasant == 239 | pleasant == 223 | pleasant == 207) = [];
% neutral(neutral == 255 | neutral == 239 | neutral == 223 | neutral == 207) = [];
% unpleasant(unpleasant == 255 | unpleasant == 239 | unpleasant == 223 | unpleasant == 207) = [];

% % Plot overal histo
% figure('color','w'); 
% subplot(3,1,1); 
% histogram(pleasant); hold on; ylabel('Pleasant','fontsize',12,'fontweight','bold')
% title(sprintf('%g preceding trials', ntrials));
% subplot(3,1,2); 
% histogram(neutral); hold on; ylabel('Neutral','fontsize',12,'fontweight','bold')
% subplot(3,1,3); 
% histogram(unpleasant); hold on; ylabel('Unpleasant','fontsize',12,'fontweight','bold')

% % Anova
% catNames = {};
% catNames(1:length(pleasant),1) = {'pleasant'};
% catNames(end+1:end+length(neutral),1) = {'neutral'};
% catNames(end+1:end+length(unpleasant),1) = {'unpleasant'};
% x = [pleasant neutral unpleasant];
% p = anova1(x,catNames);

% % Stats per subject
% % [h, ~, ~, adj_p] = fdr_bh(p_random,.05,'pdep','yes');
% h = p_random < 0.05;
% if sum(h)>0, warning(sprintf('%g subjects have a sequence that is not random for block 3!', sum(h))); end
% 
% for iSub = 1:78
%     cond1 = squeeze(Pleasant(iSub,:,:));
%     cond2 = squeeze(Neutral(iSub,:,:));
%     cond3 = squeeze(Unpleasant(iSub,:,:));
%     cond1(cond1(:,1) == 0,:) = [];
%     cond2(cond2(:,1) == 0,:) = [];
%     cond3(cond3(:,1) == 0,:) = [];
%     
%     figure('color','w');
%     subplot(3,1,1)
%     histogram(cond1); hold on; ylabel('Pleasant','fontsize',12,'fontweight','bold')
%     xlim([1 9])
%     title(sprintf('Subject %g', iSub),'fontsize',12,'fontweight','bold');
%     subplot(3,1,2)
%     histogram(cond2); hold on; ylabel('Neutral','fontsize',12,'fontweight','bold')
%     xlim([1 9])
%     subplot(3,1,3)
%     histogram(cond3); hold on; ylabel('Unpleasant','fontsize',12,'fontweight','bold')
%     xlim([1 9])
% 
%     catNames = {};
%     catNames(1:length(cond1),1) = {'pleasant'};
%     catNames(end+1:end+length(cond2),1) = {'neutral'};
%     catNames(end+1:end+length(cond3),1) = {'unpleasant'};
%     x = [cond1 cond2 cond3];
%     p = anova1(x,catNames);
% 
% 
% end



