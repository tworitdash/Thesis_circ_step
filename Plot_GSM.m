% close all;
%<<<<<<< HEAD
% c_ = load('TE_TE_Spp.mat');
c_a = load('Spp_analytical_conv.mat');

%=======
c_ = load('TE_TE_Spp.mat');
%c_a = load('TE_TE_Spp_analytical_2.mat');
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183
%c_a = load('TE_TE_S11_p_n.mat');
% GSM_ = c_.Srp;
GSM_a = c_a.Spp;

% F = 1e9:1e9:100e9;
% F_ = 4e9:0.5e9:30e9;
% F = 4e9:0.5e9:21e9;
%F = 4e9:1e9:10e9;


F = 200e9:1e9:300e9;

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
plot(F * 1e-9, db(abs(GSM_a(:, i, i)))/2, 'LineWidth', 2); grid on;
hold on;
% plot(F * 1e-9, (imag(GSM_a(:, i, i))), 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('S in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');



end


% 
% legend({'TE_{11} Numerical', 'TE_{11} Analytical', 'TE_{12} Numerical', 'TE_{12} Analytical', ...
%     'TE_{13} Numerical', 'TE_{13} Analytical',}, 'FontSize', 8, 'FontWeight', 'bold', 'location', 'southeast');
% xlim([4 15]);

%% Save 

% print('S_rp_TE', '-depsc');


