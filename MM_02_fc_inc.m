% close all;
clear;
c0 = 3e8;
% F = 20e9;

% F = 21e9:0.5e9:40e9;
% F = 4e9:0.5e9:21e9;
% F = 8e9;
% F = 10e11;
% F = 90:0.5e9:100e9;

F = 50e9;

Np = 100;
Nr = 100;

%% Modular inner cross product between the two wavegudies

X_ = load('TE_TE_Inner_P_analytical_V2.mat');
X_x = X_.X_til;
X_til = zeros(Nr, Np);

for p = 1:Np
    for r = 1:Nr
        X_til(r, p) = X_x(r, p);
    end
end



Spp = zeros(size(F, 2), Np, Np);
Spr = zeros(size(F, 2), Np, Nr);
Srp = zeros(size(F, 2), Nr, Np);
Srr = zeros(size(F, 2), Nr, Nr);


% S11pr = 

for k =  1:length(F)
    
    disp('Iteration Number:');
    disp(k);

%% Wavwguide p

rp = 0.0405319403216/2; % radius of the waveguide


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6; % Free space Permeability
erp = 1; % relative  permittivity
murp = 1; % relative Permeability
epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium


drho = rp/100;
dphi = pi/180;

[rho_, phi_] = meshgrid(eps:drho:rp, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zp = 0; 

[Qp, Zp, Yp, Kp] = QZcalculation_v2(Np, F(k), rp, erp, murp, rho_, phi_, zp, drho, dphi);

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

rr = 0.0405319403216/2.1; % radius of the waveguide
% rr = 0.0405319403216/4; % radius of the waveguide
err = 1; % relative  permittivity
murr = 1; % relative Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;
drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

[Qr, Zr, Yr, Kr] = QZcalculation_v2(Nr, F(k), rr, err, murr, rhor_, phir_, zr, drho, dphi);


Ip = eye(Np, Np);
Ir = eye(Nr, Nr);


% Qp = eye(Np(end), Np(end));
% Qr = eye(Nr(end), Nr(end));
    
X = sqrt(Kr * Zr) * X_til * sqrt(Yp * Kp); % modular inner cross product. Takes the dimension of Np \times Nr

% X = sqrt(Zr) * X_til * sqrt(Yp);

F_ = 2 * inv(Qr + X * inv(Qp) * X.');

Spp(k, :, :) = inv(Qp) * X.' * F_ * X - Ip;
Spr(k, :, :) = inv(Qp) * X.' * F_ * Qr;
Srp(k, :, :) = F_ * X;
Srr(k, :, :) = F_ * Qr - Ir;



end
% 
% save('Spp_analytical', 'Spp');
% save('Spr_analytical', 'Spr');
% save('Srp_analytical', 'Srp');
% save('Srr_analytical', 'Srr');
% % 
save('Spp_analytical_conv', 'Spp');
save('Spr_analytical_conv', 'Spr');
save('Srp_analytical_conv', 'Srp');
save('Srr_analytical_conv', 'Srr');




c_a = load('Srp_analytical_conv.mat');

GSM_a = c_a.Srp;

figure;

plot(1:1:Np, db(abs(diag(squeeze(GSM_a(1, :, :)))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Srp in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');


figure;

c_a = load('Spr_analytical_conv.mat');
GSM_a = c_a.Spr;

plot(1:1:Np, db(abs(diag(squeeze(GSM_a(1, :, :)))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spr in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

figure;

c_a = load('Spp_analytical_conv.mat');
GSM_a = c_a.Spp;

plot(1:1:Np, db(abs(diag(squeeze(GSM_a(1, :, :)))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Spp in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');

figure;

c_a = load('Srr_analytical_conv.mat');
GSM_a = c_a.Srr;

plot(1:1:Nr, db(abs(diag(squeeze(GSM_a(1, :, :)))))/2, 'LineWidth', 2); grid on;
%>>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183

xlabel('N in mode', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Srr in  dB', 'FontSize', 12, 'FontWeight', 'bold');
title(['S Parameter'], 'FontSize', 12, 'FontWeight', 'bold');
