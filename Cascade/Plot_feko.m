%%

% c_5 = load('Stt3_ratio_1_modes_5_V3.mat');
c_1 = load('Str4_ratio_1_modes_10_1mm.mat');
% c_10 = load('Srr3_ratio_1_modes_10.mat');
% c_15 = load('Srr3_ratio_1_modes_15.mat');
% c_20 = load('Srr3_ratio_1_modes_20.mat');
% c_30 = load('Stt3_ratio_1_modes_30.mat');
% c_40 = load('Srr3_ratio_1_modes_40.mat');


GSM_1 = c_1.STR;
% GSM_5 = c_5.STT;
% GSM_10 = c_10.SRR;
% GSM_15 = c_15.SRR;
% GSM_20 = c_20.SRR;
% GSM_30 = c_30.STT;
% GSM_40 = c_40.SRR;


data5 = read(rfdata.data,'4wg_touchstone_1modes_1mm.s2p');
s_params_5 = extract(data5,'S_PARAMETERS');

F = 4e9:0.5e9:21e9; % Frequency of operation
F1 = 4e9:0.5e9:35e9; % Frequency of operation Feko

figure;
plot(F * 1e-9, db(abs(squeeze(s_params_5(2, 1, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_1(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;



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

% xlim([4 35]);

% Phase_MM = atan(imag(squeeze(GSM_5(:, 1, 1)))./real(squeeze(GSM_5(:, 1, 1)))) * 180/pi;
% Phase_Feko = atan(imag(squeeze(s_params_5(1, 1, :)))./real(squeeze(s_params_5(1, 1, :)))) * 180/pi;

% Phase_MM = unwrap(angle(squeeze(GSM_5(:, 1, 1)))) * 180/pi;
% Phase_Feko = - unwrap(angle(squeeze(s_params_5(1, 1, :)))) * 180/pi;
% % 
% figure;
% plot(F1 * 1e-9, Phase_MM, 'LineWidth', 2); grid on;
% hold on;
% plot(F1 * 1e-9, Phase_Feko, 'LineWidth', 2); grid on;
% xlim([4 21]);
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
