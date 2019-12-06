close all;
c_ = load('TM_TM_Spr.mat');
c_a = load('TM_TM_Spr_analytical.mat');

GSM_ = c_.Spr;
GSM_a = c_a.Spr;

% F = 1e9:1e9:100e9;
% F = 4e9:0.5e9:30e9;
F = 4e9:0.5e9:40e9;

% F = 21e9:0.5e9:40e9;

figure
for i = 1:3


hold on;
plot(F * 1e-9, db(GSM_(:, i, i))/10, '*', 'LineWidth', 1); grid on;
hold on;
plot(F * 1e-9, db(GSM_a(:, i, i))/10, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');



end


legend({'TM_{11} Numerical', 'TM_{11} Analytical', 'TM_{12} Numerical', 'TM_{12} Analytical', ...
    'TM_{13} Numerical', 'TM_{13} Analytical',}, 'FontSize', 8, 'FontWeight', 'bold', 'location', 'northeast');
% xlim([4 30]);

print('S_pr_TM', '-depsc');