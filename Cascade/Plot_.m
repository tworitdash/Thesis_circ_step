c_5 = load('Str3_ratio_1_modes_5.mat');
c_10 = load('Str3_ratio_1_modes_10.mat');
c_15 = load('Str3_ratio_1_modes_15.mat');
c_20 = load('Str3_ratio_1_modes_20.mat');
c_30 = load('Str3_ratio_1_modes_30.mat');
c_40 = load('Str3_ratio_1_modes_40.mat');

GSM_5 = c_5.STR;
GSM_10 = c_10.STR;
GSM_15 = c_15.STR;
GSM_20 = c_20.STR;
GSM_30 = c_30.STR;
GSM_40 = c_30.STR;


figure;
plot(F * 1e-9, db(abs(squeeze(GSM_5(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_10(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_15(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_20(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_30(:, 1, 1))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(GSM_40(:, 1, 1))))/2, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('STT in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

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
