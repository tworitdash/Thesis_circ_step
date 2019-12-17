%% Inner Product Calculation
clear;


rp = 0.0405319403216/2; % radius of the waveguide P
rr = 0.0405319403216/2.1; % radius of the waveguide R
% rr = 0.0405319403216/4; % radius of the waveguide R

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


Np = 20;
Nr = 20;

Str = load('Xmn.mat');
Xmn = Str.Xmn;


drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

X_til = zeros(Nr, Np);


   for p = 1:Np
            for r = 1:Nr
                
                beta_rhop = (Xmn(p).xmn)/rp;
                beta_rhor = (Xmn(r).xmn)/rr;
                
                
                if Xmn(p).mode == "TE"
                    Nup = ((epsilonp * pi/2 .* ((Xmn(p).xmn).^2 - Xmn(p).m).^2) .* (besselj(Xmn(p).m, Xmn(p).xmn).^2)).^(-1);
                elseif Xmn(p).mode  == "TM"
                    Nup = (epsilonp .* pi/2 .* (Xmn(p).xmn).^2 .* (besselj_der(Xmn(p).m, Xmn(p).xmn)).^2).^(-1);
                end
                
                if Xmn(r).mode  == "TE"
                    Nur = ((epsilonp * pi/2 .* ((Xmn(r).xmn).^2 - Xmn(r).m).^2) .* (besselj(Xmn(p).m, Xmn(p).xmn).^2)).^(-1);
                elseif Xmn(r).mode  == "TM"
                    Nur = (epsilonp .* pi/2 .* (Xmn(r).xmn).^2 .* (besselj_der(Xmn(r).m, Xmn(r).xmn)).^2).^(-1);
                end
                
                grad_Phi_rhop = sqrt(Nup) .* cos(Xmn(p).m .* phir_) .* besselj_der(Xmn(p).m, beta_rhop .* rhor_) .* beta_rhop;
                grad_Phi_phip = (-1./rhor_) .* sqrt(Nup) .* Xmn(p).m .* sin(Xmn(p).m .* phir_) .* besselj(Xmn(p).m,  beta_rhop .* rhor_);

                grad_Phi_rhor = sqrt(Nur) .* cos(Xmn(r).m .* phir_) .* besselj_der(Xmn(r).m, beta_rhor.* rhor_) .* beta_rhor;
                grad_Phi_phir = (-1./rhor_) .* sqrt(Nur) .* Xmn(r).m .* sin(Xmn(r).m .* phir_) .* besselj(Xmn(r).m,  beta_rhor .* rhor_);

                if (Xmn(p).mode  == "TE" && Xmn(r).mode  == "TE") || (Xmn(p).mode  == "TM" && Xmn(r).mode  == "TM")
                    X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
                        .* rhor_ .* drho .* dphi;
                    X_til(r, p) = sum(sum(X_til_pr));
                elseif (Xmn(p).mode  == "TE" && Xmn(r).mode  == "TM")
                    X_til(r, p) = 0;
                elseif (Xmn(p).mode  == "TM" && Xmn(r).mode == "TE")
                    X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_rhop)...
                        .* rhor_ .* drho .* dphi;
                    X_til(r, p) = sum(sum(X_til_pr));
       
                end

%                 if (modep == "TE" && moder == "TE") || (modep == "TM" && moder == "TM")
% 
%                       A = Lommel(0, rr, beta_rhop(pm, pn), beta_rhor(rm, rn), pm - 1, rm - 1);
%                       
%                       D = Lommel(0, rr, beta_rhop(pm, pn), beta_rhor(rm, rn), pm + 1, rm + 1);
% 
%                     
%                       Icos = intphicos(0, 2*pi, pm, rm);
%                       Isin = intphisin(0, 2*pi, pm, rm);
%                       K = beta_rhop(pm, pn) .* beta_rhor(rm, rn)./4;
%                       
%                       X_til_pr =  sqrt(Nup) .* sqrt(Nur) ...
%                         .* K .* (A + D) .* (Icos + Isin);
% 
%                       X_til(pm, pn, rm, rn) = X_til_pr;                  
%                 elseif (modep == "TE" && moder == "TM")
%                     X_til(pm, pn, rm, rn) = 0;
%                elseif (modep == "TM" && moder == "TE")
%                    X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_rhop)...
%                        .* rhor_ .* drho .* dphi;
%                    X_til(pm, pn, rm, rn) = sum(sum(X_til_pr));
%        
%                 end
            end
    end
    

% csvwrite('TM_TE_Inner_P', X_til); 
save('TE_TE_Inner_P_analytical_V2', 'X_til');