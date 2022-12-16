%get HDI values from LIMO course plot difference

limo_eeg    % View results > cluster correction > course plot > select .mat result file > select EEG channel
% 58 = POz; 63 = Oz

hdi_vals = [];
index = mask(63,:) == true; %index of significant values after correction
vals = Plotted_data(index,:)';
hdi_vals(1,:) = vals(3,:);               % highest HDI
hdi_vals(2,:) = vals(2,:);               % central tendency estimate;
hdi_vals(3,:) = vals(1,:);               % lowest HDI
% hdi_vals(4,:) = vals(3,:) - vals(1,:);   % calculate difference between highest and lowest bounds

hdi_vals(4,1) = mean(hdi_vals(3,:));    %average low HDI
hdi_vals(4,2) = mean(hdi_vals(1,:));    %average high HDI

writematrix(hdi_vals, 'hdi_values.xlsx');