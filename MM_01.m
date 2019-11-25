close all;
clear;

%% Wavwguide p

mp = 1; % first digit of the mode number
Np = 1:1:3; % second digit of the mode number. p subscript is for waveguide P

modep = "TE"; % Waveguide mode polarization
% mode = "TM"

F = 20e9;

rp = 0.0405319403216/2; % radius of the waveguide


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6; % Free space Permeability
erp = 1; % relative  permittivity
murp = 1; % relative Permeability
epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium


drho = rp/1000;
dphi = pi/2000;

[rho_, phi_] = meshgrid(eps:drho:rp, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zp = 0; 

[Qp, Zp, xmn_] = QZcalculation(mp, Np, modep, F, rp, erp, murp, rho_, phi_, zp, drho, dphi);

beta_rhop = xmn_./rp;

if modep == "TE"
    Nup = (epsilonp * pi/2 .* (xmn_.^2 - mp.^2) .* (besselj(mp, xmn_)).^2).^(-1);
elseif modep == "TM"
    Nup = (epsilonp .* pi/2 .* xmn_.^2 .* besselj_der(mp, xmn_));
end



%% Wavwguide r

mr = 1; % first digit of the mode number
Nr = 1:1:5; % second digit of the mode number. p subscript is for waveguide P

moder = "TE"; % Waveguide mode polarization
% mode = "TM"

F = 1.4132e+11;

rr = 0.0405319403216/4; % radius of the waveguide
err = 1; % relative  permittivity
murr = 1; % relative Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;

drho = rr/1000;
dphi = pi/2000;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

[Qr, Zr, xmn_] = QZcalculation(mr, Nr, moder, F, rr, err, murr, rhor_, phir_, zr, drho, dphi);

beta_rhor = xmn_./rr;

if moder == "TE"
    Nur = (epsilonr * pi/2 .* (xmn_.^2 - mr.^2) .* (besselj(mr, xmn_)).^2).^(-1);
elseif moder == "TM"
    Nur = (epsilonr .* pi/2 .* xmn_.^2 .* besselj_der(mr, xmn_)).^(-1);
end

%% Modular inner cross product between the two wavegudies
X_til = zeros(Np(end), Nr(end));

for p = 1:length(Np)
    for r = 1:length(Nr)
        disp(p);
        disp(r);
        grad_Phi_rhop = Nup(p) .* cos(mp .* phir_) .* besselj_der(mp, beta_rhop(p) .* rhor_) .* beta_rhop(p);
        grad_Phi_phip = (-1./rhor_) .* Nup(p) .* mp .* sin(mp .* phir_) .* besselj(mp,  beta_rhop(p) .* rhor_);
        
        grad_Phi_rhor = Nur(r) .* cos(mr .* phir_) .* besselj_der(mr, beta_rhor(r).* rhor_) .* beta_rhor(r);
        grad_Phi_phir = (-1./rhor_) .* Nur(r) .* mr .* sin(mr .* rhor_) .* besselj(mp,  beta_rhor(r) .* rhor_);
        
        if (modep == "TE" && moder == "TE") || (modep == "TM" && moder == "TM")
            X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
                .* rhor_ .* drho .* dphi;
            X_til(p, r) = sum(sum(X_til_pr));
        elseif (modep == "TE" && moder == "TM")
            X_til(p, r) = 0;
        elseif (modep == "TM" && moder == "TE")
            X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_rhop)...
                .* rhor_ .* drho .* dphi;
            X_til(p, r) = sum(sum(X_til_pr));
        else
        end
    end
end
    
X = (Qr * Zr).^0.5 * X_til.' * (Zp\Qp).^0.5; % modular inner cross product. Takes the dimension of Np \times Nr

F_ = inv(2 * (Qr + X * pinv(Qp) * X.'));

Spp = inv(Qp) * X.' * F_ * X - eye(Np(end), Np(end));
Sws = inv(Qp) * X.' * F_ * Qr;
Ssw = F_ * X;
Sss = F_ * Qr - eye(Nr(end), Nr(end));
