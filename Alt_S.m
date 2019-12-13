close all;
clear;

% F = 20e9;
%% Mode and frequency parameters
% F = 21e9:0.5e9:40e9;
F = 4e9:0.5e9:20e9;
% F = 90:0.5e9:100e9;
mp = 1; % first digit of the mode number
Np = 1:1:5; % second digit of the mode number. p subscript is for waveguide P

mr = 1; % first digit of the mode number
Nr = 1:1:3; % second digit of the mode number. p subscript is for waveguide P


Spp = zeros(size(F, 2), Np(end), Np(end));
Spr = zeros(size(F, 2), Np(end), Nr(end));
Srp = zeros(size(F, 2), Nr(end), Np(end));
Srr = zeros(size(F, 2), Nr(end), Nr(end));

%% Waveguide parameters
    modep = "TE"; % Waveguide mode polarization
%     modep = "TM";

    rp = 0.0405319403216/2; % radius of the waveguide


    er0 = 8.85418782e-12; % Free space permittivity
    mu0 = 1.25663706e-6; % Free space Permeability
    erp = 1; % relative  permittivity
    murp = 1; % relative Permeability
    epsilonp = erp * er0;   % Permittivity in the medium
    mup = mu0 * murp;       % Permeability in the medium
    
    drho = rp/100;
    dphi = pi/200;

    [rho_, phi_] = meshgrid(eps:drho:rp, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
    zp = 0; 
    
    moder = "TE"; % Waveguide mode polarization
%     moder = "TM";

    %F = 1.4132e+11;

    rr = 0.0405319403216/2.1; % radius of the waveguide
    err = 1; % relative  permittivity
    murr = 1; % relative Permeability
    epsilonr = err * er0;   % Permittivity in the medium
    mur = mu0 * murr;

    drho = rr/100;
    dphi = pi/200;

    [rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
    zr = 0; 
   

    
for k =  1:length(F)
    
    %% Normalization constant for both P and R waveguide
    
    [Qp, Zp, Yp, xmn_p] = QZcalculation(mp, Np, modep, F(k), rp, erp, murp, rho_, phi_, zp, drho, dphi); 
    
    [Qr, Zr, Yr, xmn_r] = QZcalculation(mr, Nr, moder, F(k), rr, err, murr, rhor_, phir_, zr, drho, dphi);
    X = zeros(Nr(end), Np(end));
    Spppr = zeros(size(F, 2), Np(end), Np(end) + Nr(end));
    Srprr = zeros(size(F, 2), Nr(end), Np(end) + Nr(end));

    S = zeros(size(F, 2), Np(end) + Nr(end), Np(end) + Nr(end));
    
    for p = 1:length(Np)
        for r = 1:length(Nr)
            [Erhop, Ephip, Ezp, Hrhop, Hphip, Hzp, beta_zp] = E_and_H(rhor_, phir_, erp, murp, zp, rp, mp, Np(p), modep, F(k));
            [Erhor, Ephir, Ezr, Hrhor, Hphir, Hzr, beta_zr] = E_and_H(rhor_, phir_, err, murr, zr, rr, mr, Nr(r), moder, F(k));
            
            Xij = (Erhor .* Hphip - Hrhop .* Ephir) .* rhor_ .* drho .* dphi;
            X(r, p) = sum(sum(Xij));
        end
    end
    
Ip = eye(Np(end), Np(end));
Ir = eye(Nr(end), Nr(end));

F_ = 2 * inv(Qr + X * inv(Qp) * X.');

Spp(k, :, :) = inv(Qp) * X.' * F_ * X - Ip;
Spr(k, :, :) = inv(Qp) * X.' * F_ * Qr;
Srp(k, :, :) = F_ * X;
Srr(k, :, :) = F_ * Qr - Ir;

Spppr(k, :, :) = cat(2, squeeze(Spp(k, :, :)), squeeze(Spr(k, :, :)));
Srprr(k, :, :) = cat(2, squeeze(Srp(k, :, :)), squeeze(Srr(k, :, :)));
S(k, :, :) = cat(1, squeeze(Spppr(k, :, :)), squeeze(Srprr(k, :, :)));

squeeze(S(k, :, :))*(squeeze(S(k, :, :)))
    

    
end

save('TE_TE_Spp_analytical_2', 'Spp');
save('TE_TE_Spr_analytical_2', 'Spr');
save('TE_TE_Srp_analytical_2', 'Srp');
save('TE_TE_Srr_analytical_2', 'Srr');
