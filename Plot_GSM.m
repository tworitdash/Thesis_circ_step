close all;
c_ = load('TE_TE_Srr.mat');
c_a = load('TE_TE_Srr_analytical.mat');

GSM_ = c_.Srr;
GSM_a = c_a.Srr;

% F = 1e9:1e9:100e9;
F = 4e9:0.5e9:30e9;
F_a = 4e9:0.5e9:40e9;

% F = 21e9:0.5e9:40e9;

figure
for i = 1:3


hold on;
plot(F * 1e-9, db(GSM_(:, i, i))/10, '*', 'LineWidth', 1); grid on;
hold on;
plot(F_a * 1e-9, db(GSM_a(:, i, i))/10, 'LineWidth', 2); grid on;

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');



end


legend({'TE_{11} Numerical', 'TE_{11} Analytical', 'TE_{12} Numerical', 'TE_{12} Analytical', ...
    'TE_{13} Numerical', 'TE_{13} Analytical',}, 'FontSize', 8, 'FontWeight', 'bold', 'location', 'northeast');
xlim([4 30]);

print('S_rr_TE', '-depsc');