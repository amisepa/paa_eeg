%% Get Cohen's d effect size
clear; close all; clc;

cd 'D:\presentiment_eeg\data_clean\derivatives\LIMO_presentiment_eeg\pleasant_neutral'
tmp = load('LIMO.mat');
gpLIMO = tmp.LIMO;
eeglab; close;

%% using Betas

t = gpLIMO.data.timevect;
peakTime = find(t == -148);

for iSub = 1:78
    subPath = gpLIMO.data.data_dir{1,iSub};
    load(fullfile(subPath, 'Betas.mat')) % Beta parameter estimates

    % Difference at peak time (-148 ms) and electrode (64) using one-sample
    % ttest approach - to calculate Cohen d:
    diff(iSub,:) = Betas(64,peakTime,1) - Betas(64,peakTime,2);
    diff2(iSub,:) = Betas(64,peakTime,3) - Betas(64,peakTime,2);

end

figure('color','w')
histogram(diff)
hold on;
histogram(diff2)
legend('pleasant', 'unpleasant')

cohenD = mean(diff) / std(diff)
% cohenD = mean(diff2) / std(diff2)

%% Cohen's d on dependent groups (using HDI and t-values)

idx = find(stat_values(64,:) > 5);  % find time index of the peak effect
lowerHDI = Plotted_data(idx,1);     % lower bound of HDI of peak effect at channel 64 (O2)
upperHDI = Plotted_data(idx,3);     % higher bound of HDI of peak effect at channel 64 (O2)

% sd_diff = sqrt(N)x(upper HDI - lower HDI) / t-value
sd_diff = ( sqrt(78) * (upperHDI - lowerHDI) ) / stat_values(64,idx);
mean_diff = Plotted_data(idx,2);
cohenD = mean_diff / sd_diff;

disp(['95% HDI = [' num2str(lowerHDI) ', ' num2str(upperHDI) ']'])
disp(['Cohen''s d = ' num2str(cohenD)])


