c_ = load('TM_TM_Spr.mat');

GSM_ = c_.Spr;

% F = 1e9:1e9:100e9;
F = 1e9:0.5e9:20e9;

% F = 21e9:0.5e9:40e9;

figure;

plot(F * 1e-9, db(GSM_(:, 1, 1))/10, 'LineWidth', 2); grid on;
% plot(F * 1e-9, (GSM_(:, 1, 1)), 'LineWidth', 2); grid on;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

%legend({'Re(Z)', 'Im(Z)'}, 'FontSize', 12, 'FontWeight', 'bold');

% print('S_pp_TE_11', '-depsc');