c_5 = load('Srt3_ratio_1_modes_5.mat');
% c_10 = load('Srr3_ratio_1_modes_10.mat');
% c_15 = load('Srr3_ratio_1_modes_15.mat');
% c_20 = load('Srr3_ratio_1_modes_20.mat');
% c_30 = load('Stt3_ratio_1_modes_30.mat');
% c_40 = load('Srr3_ratio_1_modes_40.mat');

GSM_5 = c_5.SRT;
% GSM_10 = c_10.SRR;
% GSM_15 = c_15.SRR;
% GSM_20 = c_20.SRR;
% GSM_30 = c_30.STT;
% GSM_40 = c_40.SRR;


data5 = read(rfdata.data,'3wg_touchstone_5modes.s10p');
s_params_5 = extract(data5,'S_PARAMETERS');

F = 4e9:0.5e9:50e9; % Frequency of operation
F1 = 4e9:0.5e9:21e9; % Frequency of operation Feko

figure;
plot(F1 * 1e-9, db(abs(squeeze(s_params_5(1, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_5(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;



% hold on;
% plot(F * 1e-9, db(abs(squeeze(GSM_10(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(GSM_15(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(GSM_20(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(GSM_30(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(squeeze(GSM_40(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{RT} of TE_{11}, 5 modes active Feko', 'S_{RT} of TE_{11}, 5 modes active MM'},...
    'FontSize', 12, 'FontWeight', 'bold');

xlim([4 21]);


% 
% figure;
% plot(F1 * 1e-9, (angle(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, (angle(squeeze(GSM_5(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% xlim([4 21]);
% % legend({'S_{TT} of TE_{11}, 5 modes active', 'S_{TT} of TE_{11}, 10 modes active',...
%     'S_{TT} of TE_{11}, 15 modes active', 'S_{TT} of TE_{11}, 20 modes active', ...
%     'S_{TT} of TE_{11}, 30 modes active','S_{TT} of TE_{11}, 40 modes active'}, 'FontSize', 12, 'FontWeight', 'bold');

% figure;
% 
% plot(F * 1e-9, db(abs(squeeze(GSM_40(:, 1, 1))))/2 - db(abs(squeeze(GSM_10(:, 1, 1))))/2, 'LineWidth', 2); grid on;                                                                                                                                                                                       
% figure;
% plot(F * 1e-9, db(abs(squeeze(STR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('STR in  dB', 'FontSize', 12, 'FontWeight', 'bold');
% title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
% plot(F * 1e-9, db(abs(squeeze(SRT(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('SRT in  dB', 'FontSize', 12, 'FontWeight', 'bold');
% title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
% plot(F * 1e-9, db(abs(squeeze(SRR(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('SRR in  dB', 'FontSize', 12, 'FontWeight', 'bold');
% title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');