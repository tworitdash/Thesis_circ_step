
%% Data From Feko

data5 = read(rfdata.data,'S_Feko_5modes_each.s10p');
s_params_5 = extract(data5,'S_PARAMETERS');

%% Data from MM from MATLAB Script Wg2.m

% <<<<<<< HEAD
C_5 = load('Srr2_ratio_1_modes_20_fc_align.mat');
% =======
% C_5 = load('Srr2_ratio_2_modes_3_fc_align.mat');
% >>>>>>> cca7b2d7ed2d1ce0e950719caf42fcc496cc3dac

GSM_5 = C_5.Srr;

%% Frequency axis

F = 2e9:0.5e9:10e9;
F1 = 4e9:0.5e9:21e9;

%% Plot Amplitude of S parameters

figure;


% <<<<<<< HEAD
% plot(F * 1e-9, db(abs(squeeze(s_params_5(4, 4, :))))/2, 'LineWidth', 2); grid on;
% =======
plot(F1 * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
% >>>>>>> cca7b2d7ed2d1ce0e950719caf42fcc496cc3dac
hold on;
plot(F1 * 1e-9, db(abs(squeeze(GSM_5(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;


xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{RR} of TE_{11}, 5 modes active Feko', 'S_{RR} of TE_{11}, 5 modes active MM'},...
   'FontSize', 12, 'FontWeight', 'bold');
xlim([5 21]);

%% Plot Phase of S parameters

% <<<<<<< HEAD
% Phase_Feko = (angle(squeeze((s_params_5(4, 4, :))))) * 180/pi;
% =======
Phase_Feko = (angle(squeeze((s_params_5(6, 6, :))))) * 180/pi;
% >>>>>>> cca7b2d7ed2d1ce0e950719caf42fcc496cc3dac
Phase_MM =  (angle(squeeze((GSM_5(:, 1, 1))))) * 180/pi;


figure;
plot(F1 * 1e-9, Phase_Feko, 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, Phase_MM, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Phase S in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{RR} of TE_{11}, 5 modes active Feko', 'S_{RR} of TE_{11}, 5 modes active MM'},...
   'FontSize', 12, 'FontWeight', 'bold');
xlim([5 21]);

