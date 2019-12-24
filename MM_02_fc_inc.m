% close all;
clear;
c0 = 3e8;
% F = 20e9;

% F = 21e9:0.5e9:40e9;
% F = 4e9:0.5e9:21e9;
% F = 8e9;
% F = 10e11;
% F = 90:0.5e9:100e9;

% F = 200e9:1e9:300e9;

% <<<<<<< HEAD
% Np = 5;
% Nr = 5;
% =======
F = 1000e9;

Np = 1:1:50;
Np_ = 20;
Nr = 20;
Nr_ = 20;


% >>>>>>> chnages at home

%% Modular inner cross product between the two wavegudies

for i = 1:1:Np_
    
Spp = zeros(size(F, 2), i, i);
Spr = zeros(size(F, 2), i, i);
Srp = zeros(size(F, 2), i, i);
Srr = zeros(size(F, 2), i, i);
    
for k =  1:length(F)
    
disp('Iteration Number:');
disp(i);
    
    
X_ = load('TE_TE_Inner_P_analytical_V2.mat');
X_x = X_.X_til;

X_til = zeros(i, i);

for p = 1:i
    for r = 1:i
        X_til(r, p) = X_x(r, p);
    end
end

%% Wavwguide p

rp = 0.0405319403216/2; % radius of the waveguide
% rp = 0.10;


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

[Qp, Zp, Yp, Kp] = QZcalculation_v2(Np(i), F(k), rp, erp, murp, rho_, phi_, zp, drho, dphi);

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

% rr = 0.05;
err = 1; % relative  permittivity
murr = 1; % relative Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;
drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

[Qr, Zr, Yr, Kr] = QZcalculation_v2(Np(i), F(k), rr, err, murr, rhor_, phir_, zr, drho, dphi);


Ip = eye(i, i);
Ir = eye(i, i);


% Qp = eye(Np(end), Np(end));
% Qr = eye(Nr(end), Nr(end));
    
X = sqrt(Kr * Zr) * X_til * sqrt(Yp * Kp); % modular inner cross product. Takes the dimension of Np \times Nr

% X = sqrt(Zr) * X_til * sqrt(Yp);

F_ = 2 * inv(Qr + X * inv(Qp) * X.');

Spp(k, :, :) = inv(Qp) * X.' * F_ * X - Ip;
Spr(k, :, :) = inv(Qp) * X.' * F_ * Qr;
Srp(k, :, :) = F_ * X;
Srr(k, :, :) = F_ * Qr - Ir;

Spp_(i) = squeeze(Spp(k, 1, 1));
Spr_(i) = squeeze(Spr(k, 1, 1));
Srp_(i) = squeeze(Srp(k, 1, 1));
Srr_(i) = squeeze(Srr(k, 1, 1));
% c_a = load('Srp_analytical_conv.mat');
% 
% GSM_a = c_a.Srp;




end
end


save('Conv_Spp', 'Spp_');
save('Conv_spr', 'Spr_');
save('Conv_Srp', 'Srp_');
save('Conv_srr', 'Srr_');

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
% % save('Spp_analytical', 'Spp');
% % save('Spr_analytical', 'Spr');
% % save('Srp_analytical', 'Srp');
% % save('Srr_analytical', 'Srr');
% % % 
% % save('Spp_analytical_conv', 'Spp');
% % save('Spr_analytical_conv', 'Spr');
% % save('Srp_analytical_conv', 'Srp');
% % save('Srr_analytical_conv', 'Srr');

%% 



