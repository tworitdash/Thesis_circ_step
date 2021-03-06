%% Inner Product Calculation
clear;

% 
rp = 0.0405319403216/2; % radius of the waveguide P
% rr = 0.0405319403216/2.1; % radius of the waveguide R
rr = 0.0405319403216/4; % radius of the waveguide R

% rp = 0.10;
% rr = 0.05;

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


Np = 400;
Nr = 400;

Str = load('Xmn.mat');
Xmn = Str.Xmn;


drho = rr/100;
dphi = pi/180;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide
zr = 0; 

X_til = zeros(Nr, Np);


   for p = 1:Np
            for r = 1:Nr
                disp('Iteration:')
                
                disp(p);
                disp(r);
                
                xmn_p = Xmn(p).xmn;
                xmn_r = Xmn(r).xmn;
                
                beta_rhop = xmn_p/rp;
                beta_rhor = xmn_r/rr;
                
                pm = Xmn(p).m;
                rm = Xmn(r).m;
                
                modep = Xmn(p).mode;
                moder = Xmn(r).mode;
                
                
                
%                 
%                 if modep == "TE"
%                     Nup = (pi/2 .* (xmn_p.^2 - pm.^2) .* (besselj(pm, xmn_p)).^2).^(-1);
%                 elseif modep  == "TM"
%                     Nup = (pi/2 .* (xmn_p).^2 .* (besselj_der(pm, xmn_p)).^2).^(-1);
%                 end
%                 
%                 if moder  == "TE"
%                     Nur = (pi/2 .* (xmn_r.^2 - rm.^2) .* (besselj(rm, xmn_r)).^2).^(-1);
%                 elseif moder  == "TM"
%                     Nur = (pi/2 .* (Xmn(r).xmn).^2 .* (besselj_der(rm, xmn_r)).^2).^(-1);
%                 end
                Nup = 1;
                Nur = 1;
%                 
%                 if Xmn(p).m == 0
%                     
%                     grad_Phi_rhop = 0;
%                     grad_Phi_phip = 0;
%                 
%                 else
                    
                    grad_Phi_rhop = sqrt(Nup) .* cos(pm .* phir_) .* besselj_der(pm, beta_rhop .* rhor_) .* beta_rhop;
                    grad_Phi_phip = (-1./rhor_) .* sqrt(Nup) .* pm .* sin(pm .* phir_) .* besselj(pm,  beta_rhop .* rhor_);
                
%                 end
                
%                 if Xmn(r).m == 0
%                     grad_Phi_rhor = 0;
%                     grad_Phi_phir = 0;
%                 else
                    
                    grad_Phi_rhor = sqrt(Nur) .* cos(rm .* phir_) .* besselj_der(rm, beta_rhor.* rhor_) .* beta_rhor;
                    grad_Phi_phir = (-1./rhor_) .* sqrt(Nur) .* rm .* sin(rm .* phir_) .* besselj(rm,  beta_rhor .* rhor_);
%                 end
                
                if (modep  == "TE" && moder  == "TE") || (modep  == "TM" && moder  == "TM")
                    X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
                        .* rhor_ .* drho .* dphi;
                    X_til(r, p) = sum(sum(X_til_pr));
                elseif (modep  == "TE" && moder  == "TM")
                    X_til(r, p) = 0;
                elseif (modep  == "TM" && moder == "TE")
                    X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_phip)...
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
save('Inner_P_analytical_V2_ratio_2', 'X_til');