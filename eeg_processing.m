clear; close all; clc;
codeFolder = 'C:\Users\Cedric\Documents\MATLAB\presentiment_eeg\code';
dataFolder = 'C:\Users\Cedric\Documents\MATLAB\presentiment_eeg\data_raw';
outputFolder = 'C:\Users\Cedric\Documents\MATLAB\presentiment_eeg\data_clean';
load(fullfile(codeFolder, 'sInfo.mat'));
load(fullfile(codeFolder, 'chanLabels.mat'));
eeglab; close;
cd(codeFolder)

pop_editoptions('option_parallel', 0);
pop_editoptions('option_single', 0); % ensure double precision

progressbar('Processing files')
for iSub = 1:length(sInfo)
    tstart = tic;

    % Import Neuroscan data
    if sInfo(iSub).group == 1
        filepath = fullfile(dataFolder, 'group-01');
    else
        filepath = fullfile(dataFolder, 'group-02');
    end
    EEG = pop_loadcnt(fullfile(filepath, [sInfo(iSub).filename '.cnt']), ...
        'dataformat','auto','keystroke','on','memmapfile','');
    
    % Check channel labels
    if sum(contains({EEG.chanlocs.labels}, {'VEO' 'HEO'})) ~= 0
        EEG = pop_select(EEG, 'nochannel',{'HEO', 'VEO'});
    end    
    idx = contains({EEG.chanlocs.labels}, chanLabels);
    if sum(idx) ~= 64
        if EEG.nbchan > 66
            EEG = pop_select(EEG, 'nochannel',{'BP1', 'BP2', 'HL1', 'HL2'});
        else
            error(['something else is wring with this file''s channels: ' num2str(iSub) ])
        end
    end
    if ~strcmp({EEG.chanlocs(4).labels}, 'AF7')
        EEG.chanlocs.labels = chanLabels;
        warning('Check Channel labels!')
    end
    idx = find(contains({EEG.chanlocs.labels}, {'M1' 'M2' 'PO2'}));
    if ~isempty(idx)
        EEG.chanlocs(idx(1)).labels = 'TP9';
        EEG.chanlocs(idx(2)).labels = 'TP10';
        if length(idx) > 2
            EEG.chanlocs(idx(3)).labels = 'PO4';
        end
    end
    
    % Add back AFz electrode (Ground --> will be interpolated)
    EEG.nbchan = 65;
    tmp = cellstr([string({EEG.chanlocs(1:5).labels}) "AFz" string({EEG.chanlocs(6:end).labels}) ])';
    for iChan = 1:length(tmp)
        EEG.chanlocs(iChan).labels = tmp{iChan};
    end
    tmp = [ EEG.data(1:5,:); zeros(1,EEG.pnts); EEG.data(6:end,:) ];
    EEG.data = tmp;
    EEG = eeg_checkset(EEG);
    
    % Import channel locations
    locPath = fileparts(which('dipfitdefs.m'));    
    EEG = pop_chanedit(EEG,'lookup',fullfile(locPath,'standard_BEM','elec','standard_1005.elc'));
%     figure; topoplot([],EEG.chanlocs, 'style', 'blank',  'electrodes', 'labelpoint', 'chaninfo', EEG.chaninfo);

    % Now that AFz exists with locations, remove it to interpolate it later
    EEG = pop_select(EEG,'nochannel',6);

    % Downsample
    EEG = pop_resample(EEG,250);
%     oriEEG = EEG;

    % Highpass filter
%     EEG = pop_eegfiltnew(EEG, 'locutoff',0.1,'minphase',1,'plotfreqz',1);
    order = firwsord('hamming',EEG.srate,0.05);     % get filter order for transition bandwidth (e.g., 0.5 Hz)
%     gp_delay_high = (order+1-1)/2/EEG.srate     %N taps = filter length = filter order + 1 (delay in s)
    EEG = pop_firws(EEG,'forder',order,'fcutoff',0.5,'ftype','highpass', ...
         'wtype','hamming','usefftfilt',1,'minphase',1,'plotfresp',0);
    % Remove DC drifts 
%     disp('Removing DC drifts...');
%     for iChan = 1:EEG.nbchan
%         ft = fft(EEG.data(iChan,:));
%         ft(1) = 0;      % zero out the DC component
%         EEG.data(iChan,:) = ifft(ft); % Inverse transform back to time domain.
%         if mean(real(EEG.data(iChan,:))) > .005 || mean(real(EEG.data(iChan,:))) < -.005
%             warning('Mean should be closer to 0, DC drift removal must have failed.')
%         end
%     end
%     figure('color','w'); grid on
%     for iChan = 1:EEG.nbchan
%         fprintf('channel %d \n',iChan)
%         [power(iChan,:), f] = get_psd(EEG.data(iChan,:),EEG.srate*60,'hamming',50,[],EEG.srate,[0.001 1],'psd');
%         hold on; plot(f, power(iChan,:)); 
%     end

    % Reference
%     EEG = pop_reref(EEG, []);   % Average reference
    
    % Remove bad channels
    oriEEG = EEG; 
    if ~isempty(sInfo(iSub).bad_channels) 
        EEG = pop_select(EEG, 'nochannel', sInfo(iSub).bad_channels);
    end
    figure; topoplot([],EEG.chanlocs,'style','blank', 'electrodes','labelpoint','chaninfo',EEG.chaninfo);
    newPath = fullfile(outputFolder, sprintf('sub-%2.2d',iSub)); mkdir(newPath)
    saveas(gcf,fullfile(newPath, [sprintf('sub-%2.2d',iSub) '_bad-channels.png'])); close(gcf)
    
    % Detect large artifacts with ASR
%     oriEEG = EEG; 
    reconstruct = false;
    cutoff = 110;          % variance threshold (2-80)
    useriemannian = false;
    m = memory;
    maxmem = round(.9 * (m.MemAvailableAllArrays/1000000),1);  % use half of available memory in MB
    disp(['Using 90% of available memory (' num2str(round(maxmem/1000,1)) ' GB)'])
    cleanEEG = clean_asr(EEG,cutoff,[],[],[],[],[],[],[],useriemannian,maxmem);
    
    if reconstruct
        EEG = cleanEEG;
    else % remove bad segments
        mask = sum(abs(EEG.data-cleanEEG.data),1) > 1e-10;
        EEG.etc.clean_sample_mask = ~mask;
        badData = reshape(find(diff([false mask false])),2,[])';
        badData(:,2) = badData(:,2)-1;
        EEG = pop_select(EEG,'nopoint',badData);
        EEG = eeg_checkset(EEG);
%         EEG = clean_windows(EEG,0.25,[-inf 7]); 
%         vis_artifacts(EEG,oriEEG); 
    end
    
    % ICA (taking data rank into account)
    dataRank = sum(eig(cov(double(EEG.data(:,:)'))) > 1E-7);
    EEG = pop_runica(EEG,'icatype','runica','extended',1,'pca',dataRank);
%     EEG = pop_runica(EEG,'icatype','picard','pca',dataRank,'maxiter',500);
    EEG = pop_iclabel(EEG,'default');
    EEG = pop_icflag(EEG,[NaN NaN;0.95 1;0.95 1;0.99 1;0.99 1;0.99 1;NaN NaN]);
    pop_selectcomps(EEG, 1:30);
%     gong
    saveas(gcf,fullfile(newPath, [sprintf('sub-%2.2d',iSub) '_bad-components.png'])); close(gcf)
    bad_ic = find(EEG.reject.gcompreject);                      % tag bad components
    sInfo(iSub).bad_ic = bad_ic;
    warning([ 'Removing ' num2str(length(bad_ic)) ' bad components'] )
    if ~isempty(bad_ic), EEG = pop_subcomp(EEG, bad_ic, 0); end    % remove them
    
    % Interpolate bad channels
    EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical');
    pop_eegplot(EEG,1,1,1);
    saveas(gcf,fullfile(newPath, [sprintf('sub-%2.2d',iSub) '_after-ICA.fig'])); close(gcf)
    
    % Apply CSD-transformation
    locPath = fileparts(which('csd_transform.m'));
    chanlocfile = fullfile(locPath, 'chanlocs_standard_BEM_1005.ced');
    EEG = csd_transform(EEG,1,chanlocfile);  % CSD-transformation using .ced loc file
    
    % Remove trials right after key pad trials to avoid biases
    idx = find(strcmpi({EEG.event.type}, '1')) + 1;
    idx(end) = [];  % delete if it is the last trial
    EEG.event(idx) = [];

    % Epoch
    EEG = pop_epoch(EEG, {'2' '4' '8'}, [-2 2],'epochinfo', 'yes'); %2 = pleasant pictures; 4 = neutral; 8 = unpleasant; 1 = checkerboard
    
    % Remove bad trials
%     oriEEG = EEG;
    sigRMS = nan(size(EEG.data,3),1);
    sigPower = nan(size(EEG.data,3),1);
%     progressbar('Scanning epochs for bad trials')
    for iEpoch = 1:size(EEG.data,3)
        power = nan(EEG.nbchan,30);
        sigRMS(iEpoch,:) = rms(rms(EEG.data(:,:,iEpoch)));
%         sigRMS(iEpoch,:) = peak2rms(peak2rms(EEG.data(:,:,iEpoch),1));
        for iChan = 1:EEG.nbchan
            [power(iChan,:), f] = get_psd(EEG.data(iChan,:,iEpoch),EEG.srate,'hamming',50,[],EEG.srate,[70 100],'psd');
        end
%         hold on; plot(f,mean(power));
        sigPower(iEpoch,:) = rms(rms(power));  % mean across channels and RMS across frequencies
%         progressbar(iEpoch/size(EEG.data,3))
    end
%     gong    
%     badAmp = find(sigRMS > 5*std(sigRMS));        % bad epochs based on amplitude
    badAmp = find(isoutlier(sigRMS,'gesd'));        % 'mean' 'median' 'quartiles' 'grubbs' 'gesd' (default)
%     badPow = find(sigPower < 5*std(sigPower));    % bad epochs based on high-frequency power
    badPow = find(isoutlier(sigPower,'mean'));      % 'mean' (more lax) 'median' 'quartiles' 'grubbs' 'gesd'

     %%%%%%%%%% LOST CODE HERE %%%%%%%%%
   badTrials = sort([badAmp badPow]);
    EEG = pop_rejepoch(EEG, badTrials, 0);
    sinfo(iSub).nTrials = size(EEG.data,3);

%     save final file
    pop_saveset(EEG,'filepath',newpath,'filename',[sprintf('sub-%2.2d',iSub) '_eeg.set']);

    % Plot ERP
    figure; pop_timtopo(EEG, [-1450 1450], NaN);
    saveas(gcf,fullfile(newPath, [sprintf('sub-%2.2d',iSub) '_ERP.fig'])); close(gcf)

    progressbar(iSub / length(sInfo))

end

gong
save(fullfile(codeFolder, 'sInfo.mat'), 'sInfo')


%% Create STUDY

clear; close all; clc;
codeFolder = 'C:\Users\Cedric\Documents\MATLAB\presentiment_eeg\code';
dataFolder = 'C:\Users\Cedric\Documents\MATLAB\presentiment_eeg\data_clean';
load(fullfile(codeFolder, 'sInfo.mat'));
load(fullfile(codeFolder, 'chanLabels.mat'));
chanLabels(35) = {'TP9'};
chanLabels(45) = {'TP10'};
eeglab; close;

% % Remove files with less than 200 trials
% idx = [sInfo.nTrials] < 200;
% warning(['Excluding ' num2str(sum(idx)) ' files from analysis because they have less than 200 trials!'])
% sInfo(idx) = [];

commands = {};
for iSub = 1:length(sInfo)
    filepath = fullfile(dataFolder, sprintf('sub-%2.2d',iSub));
    filename = [sprintf('sub-%2.2d',iSub) '_eeg.set'];
    EEG = pop_loadset('filename',filename,'filepath',filepath);
%     figure; pop_timtopo(EEG, [-1450 1450], NaN);

    % Check channel labels (fixed)
    tmpLabels = {EEG.chanlocs.labels};
    tmpLocs = EEG.chanlocs;
    tmpData = EEG.data;
    checkLabels = strcmp(tmpLabels, chanLabels);
    if sum(checkLabels) ~= 64
        if ~strcmp(tmpLabels(4), 'AF7')
            EEG.chanlocs(4) = tmpLocs(60);
            EEG.chanlocs(5:6) = tmpLocs(4:5);
            EEG.chanlocs(7) = tmpLocs(64);
            EEG.chanlocs(8:61) = tmpLocs(6:59);
            EEG.chanlocs(62:64) = tmpLocs(61:63);
            EEG.data(4,:,:) = tmpData(60,:,:);
            EEG.data(5:6,:,:) = tmpData(4:5,:,:);
            EEG.data(7,:,:) = tmpData(64,:,:);
            EEG.data(8:61,:,:) = tmpData(6:59,:,:);
            EEG.data(62:64,:,:) = tmpData(61:63,:,:);
            EEG = eeg_checkset(EEG);
            figure; pop_timtopo(EEG, [-1450 1450], NaN);
            EEG = pop_saveset(EEG, fullfile(filepath, filename));
        end
    end
    tmpLabels = {EEG.chanlocs.labels};
    checkLabels = strcmp(tmpLabels, chanLabels);
    if sum(checkLabels) ~= 64
        disp('check this file'); break
    end
    
    EEG.saved = 'no';
    EEG.subject = sprintf('sub-%2.2d',iSub);  % need to save .set file again to update study info
    pop_saveset(EEG,'filepath',filepath,'filename',filename);
    commands = [ commands(:)' 'index' iSub 'load' fullfile(filepath, filename) ];
    [STUDY, ALLEEG] = std_editset(STUDY,ALLEEG,'name','presentiment_eeg', ...
        'commands', commands,'updatedat','on','savedat','off','rmclust','off');
    [STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG); 
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

    disp('--------------------------------------------')
    disp(['Subject ' num2str(iSub) ' done.'])
    disp('--------------------------------------------')

end

% Save
[STUDY,EEG] = pop_savestudy(STUDY,EEG,'filename','presentiment_eeg.study','filepath',dataFolder);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% Design
STUDY = std_makedesign(STUDY, ALLEEG, 1, 'name','STUDY.design 1','delfiles','off', ...
    'defaultdesign','off','variable1','type','values1',{'2','4','8'},'vartype1','categorical', ...
    'subjselect',{'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-08','sub-09', ...
    'sub-10','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19','sub-20', ...
    'sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28','sub-29','sub-30','sub-31', ...
    'sub-32','sub-33','sub-34','sub-35','sub-36','sub-37','sub-38','sub-39','sub-40','sub-41','sub-42', ...
    'sub-43','sub-44','sub-45','sub-46','sub-47','sub-48','sub-49','sub-50','sub-51','sub-52','sub-53', ...
    'sub-54','sub-55','sub-56','sub-57','sub-58','sub-59','sub-60','sub-61','sub-62','sub-63','sub-64', ...
    'sub-65','sub-66','sub-67','sub-68','sub-69','sub-70','sub-71','sub-72','sub-73','sub-74','sub-75', ...
    'sub-76','sub-77','sub-78','sub-79','sub-80','sub-81'});

% Precompute ERP
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','recompute','on','erp','on');

% Save
[STUDY,EEG] = pop_savestudy(STUDY,EEG,'filename','presentiment_eeg.study','filepath',dataFolder);
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);

% To load a .daterp file
% fileData = load('-mat', 'C:\Users\Cedric\Documents\MATLAB\presentiment_eeg\data_clean\sub-01\01.daterp');




eeglab redraw







