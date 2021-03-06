%%

c_5 = load('Srr5_ratio_1_modes_5_1mm_2cm.mat');
c_10 = load('Srr5_ratio_1_modes_5678_1mm_2cm.mat');
c_15 = load('Srr5_ratio_1_modes_57911_1mm_2cm.mat');
% c_5 = load('Srt4_ratio_1_modes_5_1mm_sl_fix.mat');
% c_5 = load('Stt5_ratio_1_modes_5_1mm_2cm.mat');
% c_15 = load('Srr3_ratio_1_modes_15.mat');
% c_20 = load('Srr3_ratio_1_modes_20.mat');
% c_30 = load('Stt3_ratio_1_modes_30.mat');
% c_40 = load('Srr3_ratio_1_modes_40.mat');


% GSM_1 = c_1.SRR;
GSM_5 = c_5.SRR;
GSM_10 = c_10.SRR;
GSM_15 = c_15.SRR;
% GSM_20 = c_20.SRR;
% GSM_30 = c_30.STT;
% GSM_40 = c_40.SRR;

% 
% data5 = read(rfdata.data,'3wg_touchstone_5modes_V2_1mm.s10p');
% s_params_5 = extract(data5,'S_PARAMETERS');

% data5 = read(rfdata.data,'4wg_touchstone_5modes_1mm.s10p');
% s_params_5 = extract(data5,'S_PARAMETERS');

data5 = read(rfdata.data,'5wg_touchstone_5modes_1mm_2cm_1mm.s10p');
s_params_5 = extract(data5,'S_PARAMETERS');



F1 = 4e9:0.5e9:50e9; % Frequency of operation
% F = 4e9:0.5e9:50e9; % Frequency of operation
F = 4e9:0.5e9:21e9; % Frequency of operation Feko

figure;


plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, db(abs(squeeze(GSM_5(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, db(abs(squeeze(GSM_10(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, db(abs(squeeze(GSM_15(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;


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

legend({'S_{TT} of TE_{11}, 5 modes active Feko', 'S_{TT} of TE_{11}, 5 modes active MM'},...
   'FontSize', 12, 'FontWeight', 'bold');

xlim([4 21]);

% Phase_MM = atan(imag(squeeze(GSM_5(:, 1, 1)))./real(squeeze(GSM_5(:, 1, 1)))) * 180/pi;
% Phase_Feko = atan(imag(squeeze(s_params_5(6, 1, :)))./real(squeeze(s_params_5(1, 1, :)))) * 180/pi;
% 
Phase_Feko = (angle(squeeze((s_params_5(6, 6, :))))) * 180/pi;
Phase_MM =  (angle(squeeze((GSM_5(:, 1, 1))))) * 180/pi;
Phase_MM_2 =  (angle(squeeze((GSM_10(:, 1, 1))))) * 180/pi;
Phase_MM_3 =  (angle(squeeze((GSM_15(:, 1, 1))))) * 180/pi;

figure;
plot(F * 1e-9, Phase_Feko, 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, Phase_MM, 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, Phase_MM_2, 'LineWidth', 2); grid on;
hold on;
plot(F1 * 1e-9, Phase_MM_3, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Phase S in deg', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter Phase'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{TT} of TE_{11}, 5 modes active Feko', 'S_{TT} of TE_{11}, 5 modes active MM'},...
   'FontSize', 12, 'FontWeight', 'bold');

xlim([4 21]);
% % % legend({'S_{TT} of TE_{11}, 5 modes active', 'S_{TT} of TE_{11}, 10 modes active',...
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
