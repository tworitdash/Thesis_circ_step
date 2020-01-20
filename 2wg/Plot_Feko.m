
%% Data From Feko

data5 = read(rfdata.data,'S_Feko_5modes_each.s10p');
s_params_5 = extract(data5,'S_PARAMETERS');

%% Data from MM from MATLAB Script Wg2.m

C_5 = load('Srr2_ratio_2_modes_3_fc_align.mat');

GSM_5 = C_5.Srr;

%% Frequency axis

F = 2e9:0.5e9:21e9;

%% Plot Amplitude of S parameters

figure;


% plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_5(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;


xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{RR} of TE_{11}, 5 modes active Feko', 'S_{RR} of TE_{11}, 5 modes active MM'},...
   'FontSize', 12, 'FontWeight', 'bold');

%% Plot Phase of S parameters

Phase_Feko = (angle(squeeze((s_params_5(6, 6, :))))) * 180/pi;
Phase_MM =  (angle(squeeze((GSM_5(:, 1, 1))))) * 180/pi;


figure;
% plot(F * 1e-9, Phase_Feko, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, Phase_MM, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Phase S in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{RR} of TE_{11}, 5 modes active Feko', 'S_{RR} of TE_{11}, 5 modes active MM'},...
   'FontSize', 12, 'FontWeight', 'bold');


