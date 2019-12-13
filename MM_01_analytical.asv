close all;
clear;
c0 = 3e8;
% F = 20e9;

% F = 21e9:0.5e9:40e9;
F = 4e9:0.5e9:40e9;
% F = 90:0.5e9:100e9;
mp = 1; % first digit of the mode number
Np = 1; % second digit of the mode number. p subscript is for waveguide P

mr = 1; % first digit of the mode number
Nr = 1; % second digit of the mode number. p subscript is for waveguide P
%% Modular inner cross product between the two wavegudies
%X_ = load('TE_TE_Inner_P_analytical.mat');
X_ = load('TE_TE_Inner_P_analytical.mat');
X_x = X_.X_til;
X_til = zeros(Nr(end), Np(end));
% X_til_f = zeros(size(F, 2), Np(end), Nr(end));
for p = 1:length(Np)
    for r = 1:length(Nr)
        disp(p);
        disp(r);
        X_til(r, p) = squeeze(X_x(mp, p, mr, r));
%         grad_Phi_rhop = sqrt(Nup(p)) .* cos(mp .* phir_) .* besselj_der(mp, beta_rhop(p) .* rhor_) .* beta_rhop(p);
%         grad_Phi_phip = (-1./rhor_) .* sqrt(Nup(p))  .* mp .* sin(mp .* phir_) .* besselj(mp,  beta_rhop(p) .* rhor_);
%         
%         grad_Phi_rhor = sqrt(Nur(r)) .* cos(mr .* phir_) .* besselj_der(mr, beta_rhor(r).* rhor_) .* beta_rhor(r);
%         grad_Phi_phir = (-1./rhor_) .* sqrt(Nur(r)) .* mr .* sin(mr .* rhor_) .* besselj(mp,  beta_rhor(r) .* rhor_);
%         
%         if (modep == "TE" && moder == "TE") || (modep == "TM" && moder == "TM")
%             X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
%                 .* rhor_ .* drho .* dphi;
%             X_til(p, r) = sum(sum(X_til_pr));
% 
%             X_til_f(k, p, r) = X_til(p, r);
%         elseif (modep == "TE" && moder == "TM")
%             X_til(p, r) = 0;
%         elseif (modep == "TM" && moder == "TE")
%             X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_rhop)...
%                 .* rhor_ .* drho .* dphi;
%             X_til(p, r) = sum(sum(X_til_pr));
%        
%         end
    end
end


Spp = zeros(size(F, 2), Np(end), Np(end));
Spr = zeros(size(F, 2), Np(end), Nr(end));
Srp = zeros(size(F, 2), Nr(end), Np(end));
Srr = zeros(size(F, 2), Nr(end), Nr(end));

Spppr = zeros(size(F, 2), Np(end), Np(end) + Nr(end));
Srprr = zeros(size(F, 2), Nr(end), Np(end) + Nr(end));

S = zeros(size(F, 2), Np(end) + Nr(end), Np(end) + Nr(end));

Zr_ = zeros(size(F, 2), Nr(end), Nr(end));
Zp_ = zeros(size(F, 2), Np(end), Np(end));
Qr_ = zeros(size(F, 2), Nr(end), Nr(end));
Qp_ = zeros(size(F, 2), Np(end), Np(end));
% S11_ = zeros(size(F, 2), Np(end), Np(end));
% S12_ = zeros(size(F, 2), Np(end), Nr(end));
% S21_ = zeros(size(F, 2), Nr(end), Np(end));
% S22_ = zeros(size(F, 2), Nr(end), Nr(end));

% S11_ = zeros(size(F, 2), Nr(end), Nr(end));
% S12_ = zeros(size(F, 2), Nr(end), Np(end));
% S21_ = zeros(size(F, 2), Np(end), Nr(end));
% S22_ = zeros(size(F, 2), Np(end), Np(end));

% S11pr = 

for k =  1:length(F)
%% Wave numbers:

 omega = 2 * pi * F(k);

 lamb = c0./F(k); % wavelength
 
 beta = 2 * pi ./ lamb;
%% Wavwguide p



modep = "TE"; % Waveguide mode polarization
% modep = "TM";



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

[Qp, Zp, Yp, xmn_] = QZcalculation(mp, Np, modep, F(k), rp, erp, murp, rho_, phi_, zp, drho, dphi);

Zp_(k, :, :) = Zp;
Qp_(k, :, :) = Qp;
% Power flow
% Pp = sqrt((Zp))/(conj(sqrt((Zp)))) * abs(Qp);

beta_rhop = xmn_./rp;

beta_zp = -1j .* sqrt(-(beta.^2 - beta_rhop.^2));

if modep == "TE"

    Kp = beta_zp ./ (omega .* mup .* epsilonp^2);

elseif modep == "TM"
    
    Kp = beta_zp ./ (omega .* mup.^2 .* epsilonp);
    
end

% if modep == "TE"
%     Nup = (epsilonp .* pi/2 .* (xmn_.^2 - mp.^2) .* (besselj(mp, xmn_)).^2).^(-1);
% elseif modep == "TM"
%     Nup = (epsilonp .* pi/2 .* xmn_.^2 .* (besselj_der(mp, xmn_)).^2).^(-1);
% end



%% Wavwguide r

moder = "TE"; % Waveguide mode polarization
% moder = "TM";

%F = 1.4132e+11;

rr = 0.0405319403216/2.1; % radius of the waveguide
err = 1; % relative  permittivity
murr = 1; % relative Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;
drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

[Qr, Zr, Yr, xmn_] = QZcalculation(mr, Nr, moder, F(k), rr, err, murr, rhor_, phir_, zr, drho, dphi);

Zr_(k, :, :) = Zr;
Qr_(k, :, :) = Qr;

beta_rhor = xmn_./rr;

beta_zr = -1j .* sqrt(-(beta.^2 - beta_rhor.^2));

if moder == "TE"

    Kr = beta_zr ./ (omega .* mur .* epsilonr^2);

elseif moder == "TM"
    
    Kr = beta_zr ./ (omega .* mur.^2 .* epsilonr);
    
end

% if moder == "TE"
%     Nur = (epsilonr * pi/2 .* (xmn_.^2 - mr.^2) .* (besselj(mr, xmn_)).^2).^(-1);
% elseif moder == "TM"
%     Nur = (epsilonr .* pi/2 .* xmn_.^2 .* (besselj_der(mr, xmn_)).^2).^(-1);
% end


% Ip = eye(Np(end), Np(end));
% Ir = eye(Nr(end), Nr(end));

Ip = eye(Np(end), Np(end));
Ir = eye(Nr(end), Nr(end));


% Qp = eye(Np(end), Np(end));
% Qr = eye(Nr(end), Nr(end));
    
X = sqrt(Kr * Zr) * X_til * sqrt(Yp * Kp); % modular inner cross product. Takes the dimension of Np \times Nr

% X = sqrt(Zr) * X_til * sqrt(Yp);

F_ = 2 * inv(Qr + X * inv(Qp) * X.');

Spp(k, :, :) = inv(Qp) * X.' * F_ * X - Ip;
%<<<<<<< HEAD

Spr(k, :, :) = inv(Qp) * X.' * F_ * Qr;
Srp(k, :, :) = F_ * X;
Srr(k, :, :) = F_ * Qr - Ir;



%% Another solution
% Ip = eye(Np(end), Np(end));
% Ir = eye(Nr(end), Nr(end));
% % 
% Spr(k, :, :) = 2 * inv(Ip + inv(Qp) * X' * inv(Qr) * X) * inv(Qp) * X';
% Spp(k, :, :) = inv(Ip + inv(Qp) * X' * inv(Qr) * X) * (inv(Qp) * X' * inv(Qr) * X - Ip);
% Srp(k, :, :) = inv(Qr) * X * (Ip - squeeze(Spp(k, :, :)));
% Srr(k, :, :) = Ir - inv(Qr) * X * squeeze(Spr(k, :, :));

% Spppr(k, :, :) = cat(2, squeeze(Spp(k, :, :)), squeeze(Spr(k, :, :)));
% Srprr(k, :, :) = cat(2, squeeze(Srp(k, :, :)), squeeze(Srr(k, :, :)));
% S(k, :, :) = cat(1, squeeze(Spppr(k, :, :)), squeeze(Srprr(k, :, :)));
% 
% squeeze(S(k, :, :))*(squeeze(S(k, :, :)))
% =======

% Spr(k, :, :) = inv(Qp) * X.' * F_ * Qr;
% Srp(k, :, :) = F_ * X;
% Srr(k, :, :) = F_ * Qr - Ir;

% Spppr(k, :, :) = cat(2, squeeze(Spp(k, :, :)), squeeze(Spr(k, :, :)));
% Srprr(k, :, :) = cat(2, squeeze(Srp(k, :, :)), squeeze(Srr(k, :, :)));
% S(k, :, :) = cat(1, squeeze(Spppr(k, :, :)), squeeze(Srprr(k, :, :)));
% 
% squeeze(S(k, :, :))*(squeeze(S(k, :, :)))

%% Another solution
% Ip = eye(Np(end), Np(end));
% Ir = eye(Nr(end), Nr(end));
% 
% Spr(k, :, :) = inv(Ip + X' * X) * inv(Qp) * X';
% Spp(k, :, :) = 2 * inv(Ip + X' * X) * (X' * X - Ip);
% Srp(k, :, :) = X * (Ip - squeeze(Spp(k, :, :)));
% Srr(k, :, :) = Ir - X * squeeze(Spr(k, :, :));
% >>>>>>> 0a21b496e0b0f28d61b4a47fcc24420c9c551183


% Me = inv(Qp) * X';
% Mh = inv(Qr) * X;
% 
% S11_(k, :, :) = inv(eye(Np(end), Np(end)) + Me * Mh) * (eye(Np(end), Np(end)) - Me * Mh);
% S21_(k, :, :) = Mh * (eye(Np(end), Np(end)) + squeeze(S11_(k, :, :)));
% S12_(k, :, :) = 2*inv(eye(Np(end), Np(end)) + Me * Mh)*Me;
% S22_(k, :, :) = (Mh * squeeze(S12_(k, :, :)) - eye(Nr(end), Nr(end)));

%% Another solution
% S11_(k, :, :) = sqrt(Qr)* inv(Qr + X * X') * (Qr - X * X') * inv(sqrt(Qr));
% S12_(k, :, :) = 2 * sqrt(Qr) * (Qr + X * X') * X * inv(sqrt(Qp));
% S21_(k, :, :) = 2 * sqrt(Qp) * (Qp + X' * X) * X' * inv(sqrt(Qr));
% S22_(k, :, :) = Qp * inv(Qp + X' * X) * (Qp - X' * X) * inv(sqrt(Qp));


% S11pr(k, :, :) = inv(X.' * X + eye(Np(end), Np(end))) * (X.' * X - eye(Np(end), Np(end)));
% S12pr(k, :, :) = 2 * inv(X.' * X + eye(Np(end), Np(end))) * X.';
% S21pr(k, :, :) = X * (eye(Np(end), Np(end)) - squeeze(S11pr(k, :, :)));
% S22pr(k, :, :) = eye(Nr(end), Nr(end)) - X * squeeze(S12pr(k, :, :));



end
% 
% save('TM_TM_Spp_analytical', 'Spp');
% save('TM_TM_Spr_analytical', 'Spr');
% save('TM_TM_Srp_analytical', 'Srp');
% save('TM_TM_Srr_analytical', 'Srr');

% save('TE_TE_Spp_analytical_2', 'Spp');
% save('TE_TE_Spr_analytical_2', 'Spr');
% save('TE_TE_Srp_analytical_2', 'Srp');
% save('TE_TE_Srr_analytical_2', 'Srr');


save('TE_TE_Spp_analytical_3', 'Spp');
save('TE_TE_Spr_analytical_3', 'Spr');
save('TE_TE_Srp_analytical_3', 'Srp');
save('TE_TE_Srr_analytical_3', 'Srr');

% 
% save('TE_TE_Spp_analytical', 'Spp');
% save('TE_TE_Spr_analytical', 'Spr');
% save('TE_TE_Srp_analytical', 'Srp');
% save('TE_TE_Srr_analytical', 'Srr');


% 
% save('TE_TE_S11', 'S11_');
% save('TE_TE_S12', 'S12_');
% save('TE_TE_S21', 'S21_');
% save('TE_TE_S22', 'S22_');

% save('TE_TE_S11_p', 'S11_');
% save('TE_TE_S12_p', 'S12_');
% save('TE_TE_S21_p', 'S21_');
% save('TE_TE_S22_p', 'S22_');

% save('X_til_TM_TM', 'X_til');
% save('X_til_TE_TE', 'X_til');
%% 
% for i = 1:3
% 
% figure;
% 
% plot(F * 1e-9, real(Qp_(:, i, i)), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, imag(Qp_(:, i, i)), 'LineWidth', 2);
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Impedance Z (\Omega)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Wave Impedance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Z)', 'Im(Z)'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
% 
% plot(F * 1e-9, real(Qr_(:, i, i)), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, imag(Qr_(:, i, i)), 'LineWidth', 2);
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Impedance Z (\Omega)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Wave Impedance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Z)', 'Im(Z)'}, 'FontSize', 12, 'FontWeight', 'bold');

% end