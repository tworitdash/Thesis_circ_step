% close all;
%<<<<<<< HEAD
% c_ = load('TE_TE_Spp.mat');
c_25 = load('Spr_ratio_1.mat');

%=======
c_1 = load('Spr_ratio_1_modes_1.mat');

c_5 = load('Spr_ratio_1_modes_5.mat');

c_7 = load('Spr_ratio_1_modes_5.mat');

c_10 = load('Srr_ratio_1_modes_10.mat');

c_30 = load('Spr_ratio_1_modes_30.mat');

% c_1_with_N = load('Spp_ratio_1_modes_1_with_N.mat');

data1 = read(rfdata.data,'S_Feko_1mode_each.s2p');
s_params_1 = extract(data1,'S_PARAMETERS');

data5 = read(rfdata.data,'S_Feko_5modes_each.s10p');
s_params_5 = extract(data5,'S_PARAMETERS');



data7 = read(rfdata.data, 'S_Feko_7modes_each.s14p');
s_params_7 = extract(data7,'S_PARAMETERS');


data8 = read(rfdata.data,'S_Feko_8modes_each.s16p');
s_params_8 = extract(data8,'S_PARAMETERS');

data9 = read(rfdata.data,'S_Feko_9modes_each.s18p');
s_params_9 = extract(data9,'S_PARAMETERS');

data10 = read(rfdata.data,'S_Feko_10modes_each.s20p');
s_params_10 = extract(data10,'S_PARAMETERS');



freq = data5.Freq;

%c_a = load('TE_TE_Spp_analytical_2.mat');
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183
%c_a = load('TE_TE_S11_p_n.mat');
% GSM_ = c_.Srp;

GSM_1 = c_1.Spr;
GSM_5 = c_5.Spr;
GSM_7 = c_5.Spr;
GSM_10 = c_10.Srr;
GSM_25 = c_25.Spr;
GSM_30 = c_30.Spr;

% GSM_1_with_N = c_1_with_N.Spp;

% F = 1e9:1e9:100e9;
% F_ = 4e9:0.5e9:30e9;
F = 4e9:0.5e9:21e9;
% F = 4e9:0.1e9:21e9;


% F = 200e9:1e9:300e9;

% F = 21e9:0.5e9:40e9;

figure;
for i = 1:1


hold on;
% plot(F * 1e-9, db(abs(GSM_(:, i, i)))/2, '*', 'LineWidth', 1); grid on;
hold on;
% plot(F * 1e-9, db((GSM_a(:, i, i))/max(GSM_a(:, i, i)))/2, 'LineWidth', 2); grid on;

%<<<<<<< HEAD
% plot(F * 1e-9, (abs(GSM_a(:, i, i))), 'LineWidth', 2); grid on;
%=======
%plot(F * 1e-9, db(abs(GSM_1(:, i, i)))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, db(abs(GSM_5(:, i, i)))/2, 'LineWidth', 2); grid on;
hold on;
% plot(freq * 1e-9, db(abs(squeeze(s_params_1(1, 1, :))))/2, '-.', 'LineWidth', 2); grid on;
hold on;
% plot(F * 1e-9, db(abs(GSM_1(:, i, i)))/2, 'LineWidth', 1); grid on;
hold on;
% plot(freq * 1e-9, db(abs(squeeze(s_params_5(1, 1, :))))/2, '--', 'LineWidth', 1); grid on;
hold on;
% plot(F * 1e-9, db(abs(GSM_5(:, i, i)))/2, 'LineWidth', 2); grid on;
hold on;
% plot(freq * 1e-9, db(abs(squeeze(s_params_7(1, 1, :))))/2, '--', 'LineWidth', 1); grid on;
hold on;
% plot(F * 1e-9, db(abs(GSM_7(:, i, i)))/2, '-.', 'LineWidth', 2); grid on;

% plot(freq * 1e-9, db(abs(squeeze(s_params_8(1, 1, :))))/2, '-.', 'LineWidth', 1); grid on;
hold on;
%plot(F * 1e-9, db(abs(GSM_8(:, i, i)))/2, '-.', 'LineWidth', 2); grid on;

% plot(freq * 1e-9, db(abs(squeeze(s_params_9(1, 1, :))))/2, '*', 'LineWidth', 1); grid on;
hold on;
%plot(F * 1e-9, db(abs(GSM_9(:, i, i)))/2, '-.', 'LineWidth', 2); grid on;

plot(freq * 1e-9, db(abs(squeeze(s_params_10(11, 11, :))))/2, 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, db(abs(GSM_10(:, i, i)))/2, '-.', 'LineWidth', 2); grid on;

% hold on;
% plot(F * 1e-9, db(abs(GSM_1_with_N(:, i, i)))/2, 'LineWidth', 2); grid on;
hold on;
%plot(F * 1e-9, -db(abs(GSM_10(:, i, i))), 'LineWidth', 2); grid on;
hold on;
% plot(F * 1e-9, db(abs(GSM_30(:, i, i)))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183





end
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

% 
% legend({'|S_{pp}| of TE_{11} mode with 1 mode on each side',...
%     '|S_{pp}| of TE_{11} mode with 10 modes on each side',...
%     '|S_{pp}| of TE_{11} mode with 25 modes on each side', ...
%     '|S_{pp}| of TE_{11} mode with 30 modes on each side'},...
%     'FontSize', 12, 'FontWeight', 'bold', 'location', 'southeast');

% legend({'|S_{pp}| of TE_{11} with 1 mode active - Feko',...
%     '|S_{pp}| of TE_{11} with 1 mode active - MM', ...
%     '|S_{pp}| of TE_{11} with 5 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 5 modes active - MM',...
%     '|S_{pp}| of TE_{11} with 7 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 7 modes active - MM', ...
%     '|S_{pp}| of TE_{11} with 8 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 8 modes active - MM'...
%     '|S_{pp}| of TE_{11} with 9 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 9 modes active - MM'...
%     '|S_{pp}| of TE_{11} with 10 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 10 modes active - MM'},...
%     'FontSize', 12, 'FontWeight', 'bold', 'location', 'southeast');


% legend({'|S_{pp}| of TE_{11} with 1 mode active - Feko',...
%     '|S_{pp}| of TE_{11} with 5 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 7 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 8 modes active - Feko'...
%     '|S_{pp}| of TE_{11} with 9 modes active - Feko',...
%     '|S_{pp}| of TE_{11} with 10 modes active - Feko'},...
%     'FontSize', 12, 'FontWeight', 'bold', 'location', 'southeast');


legend({'|S_{rr}| of TE_{11} with 10 modes active - Feko',...
    '|S_{rr}| of TE_{11} with 10 modes active - MM'},...
    'FontSize', 12, 'FontWeight', 'bold', 'location', 'southeast');
% xlim([4.56 21]); 

% print('Spp_conv_ratio_1', '-depsc');

%% 



Np_ = 25;

figure(1);
c_a = load('Conv_Srp_1_mode_1.mat');
Srp_ = c_a.Srp_;
% 
% 
plot(1:1:Np_, db(abs((squeeze(Srp_))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Srp in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

print('Srp_radii_ratio_2_modes_ratio_1', '-depsc');


figure(2);

c_a = load('Conv_Spr_1_mode_1.mat');
Spr_ = c_a.Spr_;

plot(1:1:Np_, db(abs((squeeze(Spr_))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spr in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

print('Spr_radii_ratio_2_modes_ratio_1', '-depsc');

figure(3);

c_a = load('Conv_Spp_1_mode_1.mat');
Spp_ = c_a.Spp_;

plot(1:1:Np_, db(abs((squeeze(Spp_))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spp in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

print('Spp_radii_ratio_2_modes_ratio_1', '-depsc');

figure(4);

c_a = load('Conv_Srr_1_mode_1.mat');
Srr_ = c_a.Srr_;

plot(1:1:Np_, db(abs((squeeze(Srr_))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Srr in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

print('Srr_radii_ratio_2_modes_ratio_1', '-depsc');