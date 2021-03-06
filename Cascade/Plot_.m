c_5 =  load('Stt3_ratio_1_modes_5_1mm.mat');
c_10 = load('Stt3_ratio_1_modes_10_1mm.mat');
c_15 = load('Stt3_ratio_1_modes_15_1mm.mat');
c_20 = load('Stt3_ratio_1_modes_20_1mm.mat');
c_25 = load('Stt3_ratio_1_modes_25_1mm.mat');
% c_30 = load('Stt3_ratio_1_modes_30_2cm.mat');
% c_40 = load('Srr3_ratio_1_modes_40.mat');

GSM_5 = c_5.STT;
GSM_10 = c_10.STT;
GSM_15 = c_15.STT;
GSM_20 = c_20.STT;
GSM_25 = c_25.STT;
% GSM_30 = c_30.STT;
% GSM_40 = c_40.SRR;
F = 4e9:0.5e9:50e9; % Frequency of operation

% F = 4e9:0.5e9:50e9; % Frequency of operation

figure;
plot(F * 1e-9, db(abs(squeeze(GSM_5(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_10(:, 1, 1))))/2, 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_15(:, 1, 1))))/2, '-.', 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_20(:, 1, 1))))/2, 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_25(:, 1, 1))))/2, 'LineWidth', 1); grid on;
% % hold on;
% plot(F * 1e-9, db(abs(squeeze(GSM_40(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('STT in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{TT} of TE_{11}, 5 modes active', 'S_{TT} of TE_{11}, 10 modes active',...
    'S_{TT} of TE_{11}, 15 modes active', 'S_{TT} of TE_{11}, 20 modes active', ...
     'S_{TT} of TE_{11}, 25 modes active'}, 'FontSize', 8, 'FontWeight', 'bold');

 figure;
plot(F * 1e-9, angle(squeeze(GSM_5(:, 1, 1))) * 180/pi, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(GSM_10(:, 1, 1)))) * 180/pi, 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(GSM_15(:, 1, 1)))) * 180/pi, '-.', 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(GSM_20(:, 1, 1)))) * 180/pi, 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, angle((squeeze(GSM_25(:, 1, 1)))) * 180/pi, 'LineWidth', 1); grid on;
% % hold on;
% plot(F * 1e-9, db(abs(squeeze(GSM_40(:, 1, 1))))/2, '-.', 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Phase of STT [deg]', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'S_{TT} of TE_{11}, 5 modes active', 'S_{TT} of TE_{11}, 10 modes active',...
    'S_{TT} of TE_{11}, 15 modes active', 'S_{TT} of TE_{11}, 20 modes active', ...
     'S_{TT} of TE_{11}, 25 modes active'}, 'FontSize', 8, 'FontWeight', 'bold');
 
 
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
