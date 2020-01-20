function [Spp, Spr, Srp, Srr] = GSM(Nr, Np, F, rp, rr, erp, murp, err, murr, X_)

% close all;


Spp = zeros(length(Np), length(Np));
Spr = zeros(length(Np), length(Nr));
Srp = zeros(length(Nr), length(Np));
Srr = zeros(length(Nr), length(Nr));

%% Frequency independent Modular inner cross product between the two wavegudies




 X_til = zeros(length(Nr), length(Np));

 for p = 1:length(Np)
     for r = 1:length(Nr)
           X_til(r, p) = X_(r, p);
     end
 end
    

% disp('Iteration Number:');
% disp(k);
    
    

%% Wavwguide p

drho = rp/100;
dphi = pi/180;

[rho_, phi_] = meshgrid(eps:drho:rp, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zp = 0; 

[Qp, Zp, Yp, Kp] = QZcalculation_v2(Np, F, rp, erp, murp, rho_, phi_, zp, drho, dphi);

% figure;
% plot(1:1:Np, (real(diag(Qp))), 'LineWidth', 2); grid on;
% hold on;
% plot(1:1:Np, (imag(diag(Qp))), 'LineWidth', 2); grid on;
% 
% xlabel('n in TE_{m, 2} modes', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Normalization Constant Q_{m, 2}', 'FontSize', 12, 'FontWeight', 'bold');
% %title(['Normalization Constant for', mode,'_{m, 2} modes'], 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Q)', 'Im(Q)'}, 'FontSize', 12, 'FontWeight', 'bold');

% 
% 
% 
%% Wavwguide r

drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

[Qr, Zr, Yr, Kr] = QZcalculation_v2(Nr, F, rr, err, murr, rhor_, phir_, zr, drho, dphi);


Ip = eye(length(Np), length(Np));
Ir = eye(length(Nr), length(Nr));


% Qp = eye(Np(end), Np(end));
% Qr = eye(Nr(end), Nr(end));
    
X = sqrt(Kr * Zr) * X_til * sqrt(Yp * Kp); % modular inner cross product. Takes the dimension of Np \times Nr

% X = sqrt(Zr) * X_til * sqrt(Yp);

F_ = 2 * inv(Qr + X * inv(Qp) * X.');

Spp = inv(Qp) * X.' * F_ * X - Ip;
Spr = inv(Qp) * X.' * F_ * Qr;
Srp = (F_ * X)';
Srr = F_ * Qr - Ir;

% Spp_(i) = squeeze(Spp(k, 1, 1));
% Spr_(i) = squeeze(Spr(k, 1, 1));
% Srp_(i) = squeeze(Srp(k, 1, 1));
% Srr_(i) = squeeze(Srr(k, 1, 1));
% c_a = load('Srp_analytical_conv.mat');
% 
% GSM_a = c_a.Srp;

end


% save('Conv_Spp_16', 'Spp_');
% save('Conv_spr_16', 'Spr_');
% save('Conv_Srp_16', 'Srp_');
% save('Conv_srr_16', 'Srr_');

%% 

% figure(1);
% 
% plot(1:1:Np_, db(abs((squeeze(Srp_))))/2, 'LineWidth', 2); grid on;
% %>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183
% 
% xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Srp in  dB', 'FontSize', 12, 'FontWeight', 'bold');
% title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');
% 
% 
% figure(2);
% 
% % c_a = load('Spr_analytical_conv.mat');
% % GSM_a = c_a.Spr;
% 
% plot(1:1:Np_, db(abs((squeeze(Spr_))))/2, 'LineWidth', 2); grid on;
% %>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183
% 
% xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Spr in  dB', 'FontSize', 12, 'FontWeight', 'bold');
% title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure(3);
% 
% % c_a = load('Spp_analytical_conv.mat');
% % GSM_a = c_a.Spp;
% 
% plot(1:1:Np_, db(abs((squeeze(Spp_))))/2, 'LineWidth', 2); grid on;
% %>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183
% 
% xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Spp in  dB', 'FontSize', 12, 'FontWeight', 'bold');
% title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure(4);
% 
% % c_a = load('Srr_analytical_conv.mat');
% % GSM_a = c_a.Srr;
% 
% plot(1:1:Np_, db(abs((squeeze(Srr_))))/2, 'LineWidth', 2); grid on;
% %>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183
% 
% xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Srr in  dB', 'FontSize', 12, 'FontWeight', 'bold');
% title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');
% 
% % 
% save('Spp_ratio_1_modes_7', 'Spp');
% save('Spr_ratio_1_modes_7', 'Spr');
% save('Srp_ratio_1_modes_7', 'Srp');
% save('Srr_ratio_1_modes_7', 'Srr');
% % % 
% % save('Spp_analytical_conv', 'Spp');
% % save('Spr_analytical_conv', 'Spr');
% % save('Srp_analytical_conv', 'Srp');
% % save('Srr_analytical_conv', 'Srr');
%% 
% 
% figure;
% plot(F, db(abs(squeeze(Spp(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% figure;
% plot(F, db(abs(squeeze(Spr(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% figure;
% plot(F, db(abs(squeeze(Srp(:, 1, 1))))/2, 'LineWidth', 2); grid on;
% 
% figure;
% plot(F, db(abs(squeeze(Srr(:, 1, 1))))/2, 'LineWidth', 2); grid on;

%% 


