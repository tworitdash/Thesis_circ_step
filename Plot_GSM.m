close all;
c_ = load('TE_TE_Srr_analytical.mat');

GSM_ = c_.Srr;

% F = 1e9:1e9:100e9;
F = 4e9:0.5e9:30e9;

% F = 21e9:0.5e9:40e9;

figure
for i = 1:3


hold on;
plot(F * 1e-9, db(GSM_(:, i, i))/10, 'LineWidth', 2); grid on;
hold on;
% plot(F * 1e-9, (GSM_(:, 1, 1)), 'LineWidth', 2); grid on;
xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

legend({'TE_{11}', 'TE_{12}', 'TE_{13}'}, 'FontSize', 12, 'FontWeight', 'bold');

end
print('S_rr_TE', '-depsc');