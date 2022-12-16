clear;clc;close all;
dataFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\source_data\group-02';
codeFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\code';
cleanDataFolder = 'C:\Users\IONSLAB\Documents\RESEARCH\precog_v1\clean_data2';
load(fullfile(codeFolder, 'sInfo3.mat'));
cd(dataFolder)
eeglab
chanLocs = fileparts(which('dipfitdefs.m'));
iSubj = 2;
EEG = pop_loadcnt(fullfile(dataFolder, [sInfo(iSubj).filename '.cnt']), 'dataformat', 'auto', 'memmapfile', '');
EEG = pop_select(EEG, 'nochannel',{'VEO','HEO'});
EEG.chanlocs(35).labels = 'TP9';
EEG.chanlocs(45).labels = 'TP10';
EEG = pop_chanedit(EEG, 'lookup',fullfile(chanLocs, 'standard_BESA', 'standard-10-5-cap385.elp'));
EEG = eeg_checkset(EEG);
EEG = pop_resample(EEG, 250);
EEG = eeg_checkset(EEG);
oriEEG = EEG;
EEG = pop_select(EEG, 'nochannel', sInfo(iSubj).bad_channels);
EEG = eeg_checkset(EEG);
EEG = pop_interp(EEG, oriEEG.chanlocs, 'spherical');
EEG = eeg_checkset(EEG);

%raw
EEG0 = EEG;

%Filter specs
fcutoff_high = 1;              %high-pass > 0.5 Hz
fcutoff_low = 50;               %low-pass < 50 Hz
trans_highpass = 0.5;             %transition bandwidth low-pass (in Hz)
trans_lowpass = 12.5;           %transition bandwidth high-pass (in Hz)
order_highpass = firwsord('hamming', EEG.srate, trans_highpass)     %get filter order (higher the more accurate but the more group delay if using linear causal filter)
order_lowpass = firwsord('hamming', EEG.srate, trans_lowpass)       %get filter order (higher the more accurate but the more group delay if using linear causal filter)
gp_delay_high = (order_highpass+1-1)/2/EEG.srate     %N taps = filter length = filter order + 1 (delay is provided in seconds)
gp_delay_low = (order_lowpass+1-1)/2/EEG.srate     %in seconds when dividing by sample rate at the end
% gp_delay_low = (order_low+1-1)/2     %in samples
b_high = firws(order_highpass, fcutoff_high/(EEG.srate/2), 'high', windows('hamming', order_highpass + 1)); %Unit for f is pi rad/sample (i.e. normalized to Nyquist, MATLAB standard).
b_low = firws(order_lowpass, fcutoff_low/(EEG.srate/2), 'low', windows('hamming', order_lowpass + 1)); %Unit for f is pi rad/sample (i.e. normalized to Nyquist, MATLAB standard).
% figure; plot(b_high); title(['high-pass coefficients; order = ' num2str(order_highpass) '; delay = ' num2str(gp_delay_high) ]);
% figure; plot(b_low); title(['low-pass coefficients; order = ' num2str(order_lowpass) '; delay = ' num2str(gp_delay_low) ]);

delay_clever_bandpass_filter = (max([order_lowpass+1, order_highpass+1]) - 1)/2/EEG.srate

%non-causal zero-phase
EEG1 = EEG;
EEG1.data = fir_filterdcpadded(b_high, 1, EEG1.data', 0)';     %best method for linear filtering
EEG1.data = fir_filterdcpadded(b_low, 1, EEG1.data', 0)';    %best method for linear filtering

%causal zero-phase
EEG2 = EEG;
EEG2.data = fir_filterdcpadded(b_high, 1, EEG2.data', 1)';     %best method for linear filtering
EEG2.data = fir_filterdcpadded(b_low, 1, EEG2.data', 1)';    %best method for linear filtering
% EEG3.data = filter(b_high3, 1, EEG3.data); % linear causal filter
% EEG3.data = filter(b_low3, 1, EEG3.data); % linear causal filter

%Causal minimum-phase
EEG3 = EEG;
EEG3 = pop_firws(EEG3, 'forder', order_highpass, 'fcutoff', fcutoff_high, 'ftype', 'highpass', 'wtype', 'hamming', 'minphase', true, 'plotfresp', false);
EEG3 = pop_firws(EEG3, 'forder', order_lowpass, 'fcutoff', fcutoff_low, 'ftype', 'lowpass', 'wtype', 'hamming', 'minphase', true, 'plotfresp', false);

% Causal minim-phase Butterworth order 10
[b,a]= butter(10,lowcutoff/nf,'high');   % only Butterworth highpass for scale 1
signal = filtfilt(b,a,signal);

% Remove DC offset only
ft = fft(EEG.data(iChan,:));
ft(1) = 0;  %zero out the DC component
EEG.data(iChan,:) = ifft(ft); % Inverse transform back to time domain.

%Average ref
EEG0 = pop_reref(EEG0, []);
EEG1 = pop_reref(EEG1, []);
EEG2 = pop_reref(EEG2, []);
EEG3 = pop_reref(EEG3, []);

%Segment
EEG0 = pop_epoch(EEG0, {'2'}, [-1.5 1]);
EEG1 = pop_epoch(EEG1, {'2'}, [-1.5 1]);
EEG2 = pop_epoch(EEG2, {'2'}, [-1.5 1]);
EEG3 = pop_epoch(EEG3, {'2'}, [-1.5 1]); %2 = pleasant pictures; 4 = neutral; 8 = unpleasant; 

%power spectral density
[pxx0, f0] = pwelch(EEG0.data(58,:), EEG.srate, EEG.srate/2, [], EEG.srate, 'psd');
[pxx1, f1] = pwelch(EEG1.data(58,:), EEG.srate, EEG.srate/2, [], EEG.srate, 'psd');
[pxx2, f2] = pwelch(EEG2.data(58,:), EEG.srate, EEG.srate/2, [], EEG.srate, 'psd');
[pxx3, f3] = pwelch(EEG3.data(58,:), EEG.srate, EEG.srate/2, [], EEG.srate, 'psd');

%plot PSD
figure; plot(f0,10*log10(pxx0)); hold on; plot(f1,10*log10(pxx1)); hold on; plot(f2,10*log10(pxx2)); hold on; plot(f3,10*log10(pxx3));
xlim([0 80]); xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)')
legend('raw', 'zero-phase non-causal', 'zero-phase causal', 'min-phase causal');
title(['high-pass: ' num2str(fcutoff_high) '; low-pass: ' num2str(fcutoff_low)]);

%plot 1 epoch 1 channel
figure; plot(EEG0.data(58,:,10)); hold on; plot(EEG1.data(58,:,10)); hold on; plot(EEG2.data(58,:,10)); hold on; plot(EEG3.data(58,:,10));
legend('raw', 'zero-phase non-causal', 'zero-phase causal', 'min-phase causal');

%Plot differences in raw data: channel POz, trial 10
figure; plot(EEG0.data(58,:,10)); hold on; plot(EEG1.data(58,:,10)); hold on; plot(EEG0.data(58,:,10) - EEG1.data(58,:,10)); 
legend('raw', 'zero-phase non-causal', 'difference'); title(['high-pass: ' num2str(fcutoff_high) '; low-pass: ' num2str(fcutoff_low)]);
figure; plot(EEG0.data(58,:,10)); hold on; plot(EEG2.data(58,:,10)); hold on; plot(EEG0.data(58,:,10) - EEG2.data(58,:,10)); 
legend('raw', 'zero-phase causal', 'difference'); title(['high-pass: ' num2str(fcutoff_high) '; low-pass: ' num2str(fcutoff_low)]);
figure; plot(EEG0.data(58,:,10)); hold on; plot(EEG3.data(58,:,10)); hold on; plot(EEG0.data(58,:,10) - EEG3.data(58,:,10)); 
legend('raw', 'min-phase causal', 'difference'); title(['high-pass: ' num2str(fcutoff_high) '; low-pass: ' num2str(fcutoff_low)]);
figure; plot(EEG1.data(58,:,10)); hold on; plot(EEG2.data(58,:,10)); hold on; plot(EEG1.data(58,:,10) - EEG2.data(58,:,10)); 
legend('zero-phase non-causal', 'zero-phase causal', 'difference'); title(['high-pass: ' num2str(fcutoff_high) '; low-pass: ' num2str(fcutoff_low)]);
figure; plot(EEG1.data(58,:,10)); hold on; plot(EEG3.data(58,:,10)); hold on; plot(EEG1.data(58,:,10) - EEG3.data(58,:,10)); 
legend('zero-phase non-causal', 'min-phase causal', 'difference'); title(['high-pass: ' num2str(fcutoff_high) '; low-pass: ' num2str(fcutoff_low)]);
figure; plot(EEG2.data(58,:,10)); hold on; plot(EEG3.data(58,:,10)); hold on; plot(EEG2.data(58,:,10) - EEG3.data(58,:,10)); 
legend('zero-phase causal', 'min-phase causal', 'difference'); title(['high-pass: ' num2str(fcutoff_high) '; low-pass: ' num2str(fcutoff_low)]);

%ERP plots: average and topo
figure; pop_timtopo(EEG0, [-500 500], 100, 'raw');
figure; pop_timtopo(EEG1, [-500 500], 100, 'zero-phase non-causal');
figure; pop_timtopo(EEG2, [-500 500], 100, 'zero-phase causal');
figure; pop_timtopo(EEG3, [-500 500], 100, 'min-phase causal');

%Time-frequency plots
figure; pop_erpimage(EEG0, 1, 58,[],'raw',10,1,{'2'},[],'type' ,'yerplabel','\muV','erp','on','limits',[-500 500 NaN NaN NaN NaN NaN NaN] ,'cbar','on');
figure; pop_erpimage(EEG1, 1, 58,[],'zero-phase non-causal',10,1,{'2'},[],'type' ,'yerplabel','\muV','erp','on','limits',[-500 500 NaN NaN NaN NaN NaN NaN] ,'cbar','on');
figure; pop_erpimage(EEG2, 1, 58,[],'zero-phase causal',10,1,{'2'},[],'type' ,'yerplabel','\muV','erp','on','limits',[-500 500 NaN NaN NaN NaN NaN NaN] ,'cbar','on');
figure; pop_erpimage(EEG3, 1, 58,[],'min-phase causal',10,1,{'2'},[],'type' ,'yerplabel','\muV','erp','on','limits',[-500 500 NaN NaN NaN NaN NaN NaN] ,'cbar','on');

% figure; pop_erpimage(EEG, 1, [1], [[]], 'min-phase + causal low-pass 40 Hz',10,1,{},[],'' ,'yerplabel','\muV','erp','on','limits',[-1500 0 NaN NaN NaN NaN NaN NaN] ,'cbar','on','topo', { [1] EEG.chanlocs EEG.chaninfo } );
% figure; pop_timtopo(EEG, [-1000 998], [NaN], 'min-phase + causal low-pass 40 Hz');
% figure; pop_plottopo(EEG, [1:64] , 'min-phase + causal low-pass 40 Hz', 0, 'ydir',1);

%%
%group delay = (N_taps – 1) / (2 * sample_rate)     %filter length = N taps
%phase shift = (filter length - 1) / 2 samples (i.e. group delay)
%filter order = number of filter coefficients / filter length -1
%filter length = N taps = filter order + 1.

% Example:
sig = [ 0 0 0 0 0 1 0 0 0 0 0 ]; % test signal (impulse)
b_high = [ 1 1 1 1 1 ] / 5; % some crude boxcar filter for demonstration purposes only, linear-phase, length = 5, order = 4, group delay = 2
fsig = filter( b_high, 1, sig ) % causal filter
figure; plot( -5:5, [ sig; fsig ]') % the filtered impulse in the output does not start before the impulse in the input

fsig = filter( b_high, 1, [ sig 0 0 ] ) % padded causal filter
fsig = fsig( 3:end ) % delay correction by group delay, this is what makes the filter non-causal and zero-phase
figure; plot( -5:5, [ sig; fsig ]') % the filtered impulse in the output starts before the impulse in the input BUT everything before x = -2 is unaffected

b_high = minphaserceps( firws( 10, 1/5, 'high' ) ); % minimum-phase non-linear high-pass
fsig_high = filter( b_high, 1, sig ); % causal filter
figure; plot( -5:5, [ sig; fsig_high ]') % Causality preserved :)

fsig_high_low = filter( b_high, 1, [ fsig_high 0 0 ] ); % the low-pass from above
fsig_high_low = fsig_high_low( 3:13 ); % delay correction/zero-phase
figure; plot( -5:5, [ sig; fsig_high_low ]') % Causality violated :(
