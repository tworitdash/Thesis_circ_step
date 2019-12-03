
%% Inner Product Calculation
clear;
rp = 0.0405319403216/2; % radius of the waveguide P
rr = 0.0405319403216/4; % radius of the waveguide R


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

erp = 1; % relative  permittivity
murp = 1; % relative Permeability
epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium
err = 1; % relative  permittivity
murr = 1; % relative Permeability
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;

modep = "TE";
% modep = "TM";
moder = "TE";
% moder = "TM";

Mp = 1:1:5;
Np = 1:1:5;
Mr = 1:1:5;
Nr = 1:1:5;



x_TE = csvread('TE_Xmn');
x_TM = csvread('TM_Xmn');

if modep == "TE"
    xmn_p = x_TE;
elseif modep == "TM"
    xmn_p = x_TM;
end

if moder == "TE"
    xmn_r = x_TE;
elseif moder == "TM"
    xmn_r = x_TM;
end

beta_rhop = xmn_p./rp;

beta_rhor = xmn_r./rr;

drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

X_til = zeros(Mp(end), Np(end), Mr(end), Nr(end));

for pm = 1:length(Mp)
    for pn = 1:length(Np)
        for rm = 1:length(Mr)
            for rn = 1:length(Nr)
                
                if modep == "TE"
                    Nup = (epsilonp * pi/2 .* (xmn_p(pm, pn).^2 - pm.^2) .* (besselj(pm, xmn_p(pm, pn))).^2).^(-1);
                elseif modep == "TM"
                    Nup = (epsilonp .* pi/2 .* xmn_p(pm, pn).^2 .* (besselj_der(pm, xmn_p(pm, pn))).^2);
                end
                
                if moder == "TE"
                    Nur = (epsilonr * pi/2 .* (xmn_r(rm, rn).^2 - rm.^2) .* (besselj(rm, xmn_r(rm, rn))).^2).^(-1);
                elseif moder == "TM"
                    Nur = (epsilonr .* pi/2 .* xmn_r(rm, rn).^2 .* (besselj_der(rm, xmn_r(rm, rn))).^2).^(-1);
                end
                
%                 disp(p);
%                 disp(r);
%                 grad_Phi_rhop = sqrt(Nup) .* cos(pm .* phir_) .* besselj_der(pm, beta_rhop(pm, pn) .* rhor_) .* beta_rhop(pm, pn);
%                 grad_Phi_phip = (-1./rhor_) .* sqrt(Nup) .* pm .* sin(pm .* phir_) .* besselj(pm,  beta_rhop(pm, pn) .* rhor_);
% 
%                 grad_Phi_rhor = sqrt(Nur) .* cos(rm .* phir_) .* besselj_der(rm, beta_rhor(rm, rn).* rhor_) .* beta_rhor(rm, rn);
%                 grad_Phi_phir = (-1./rhor_) .* sqrt(Nur) .* rm .* sin(rm .* rhor_) .* besselj(rm,  beta_rhor(rm, rn) .* rhor_);

                if (modep == "TE" && moder == "TE") || (modep == "TM" && moder == "TM")

                      A = Lommel(0, rr, beta_rhop(pm, pn), beta_rhor(rm, rn), pm - 1, rm - 1);
                      B = Lommel(0, rr, beta_rhop(pm, pn), beta_rhor(rm, rn), pm - 1, rm + 1);
                      C = Lommel(0, rr, beta_rhop(pm, pn), beta_rhor(rm, rn), pm + 1, rm - 1);
                      D = Lommel(0, rr, beta_rhop(pm, pn), beta_rhor(rm, rn), pm + 1, rm + 1);

%                        A = Lommel(0, rr, beta_rhor(rm, rn), beta_rhop(pm, pn), pm - 1, rm - 1);
%                        B = Lommel(0, rr, beta_rhor(rm, rn), beta_rhop(pm, pn), pm - 1, rm + 1);
%                        C = Lommel(0, rr, beta_rhor(rm, rn), beta_rhop(pm, pn), pm + 1, rm - 1);
%                        D = Lommel(0, rr, beta_rhor(rm, rn), beta_rhop(pm, pn), pm + 1, rm + 1);
                      I1rho = A - B - C + D;
                      I1phi = (beta_rhop(pm, pn) .* beta_rhor(rm, rn))./(4 .* pm .* rm) .* (A + B + C + D);
                      I2rho = intphicos(0, 2*pi, pm, rm);
                      I2phi = intphisin(0, 2*pi, pm, rm);
                      
                      X_til_pr =  sqrt(Nup) .* sqrt(Nur) ...
                        .* (beta_rhop(pm, pn) .* beta_rhor(rm, rn) .* I1rho .* I2rho + pm .* rm .* I1phi .* I2phi);

                      X_til(pm, pn, rm, rn) = X_til_pr;
%                     X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
%                         .* rhor_ .* drho .* dphi;
%                     X_til(pm, pn, rm, rn) = sum(sum(X_til_pr));
                elseif (modep == "TE" && moder == "TM")
                    X_til(pm, pn, rm, rn) = 0;
                elseif (modep == "TM" && moder == "TE")
%                     X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_rhop)...
%                         .* rhor_ .* drho .* dphi;
%                     X_til(pm, pn, rm, rn) = sum(sum(X_til_pr));
       
                end
            end
        end
    end
    
 
end

% csvwrite('TM_TE_Inner_P', X_til); 
save('TE_TE_Inner_P_analytical', 'X_til');