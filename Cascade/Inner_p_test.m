function [X_til] = Inner_p_test(Nr, Np, rp, rr, erp, murp, err, murr) 

%% Inner Product Calculation


er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;


Str = load('Xmn.mat');
Xmn = Str.Xmn;


drho = rr/1000;
dphi = pi/1800;

[rhor_, phir_] = meshgrid(eps:drho:rr, eps:dphi:3/2*pi-eps);  % domain for the fields on one cross-section of the waveguide

X_til = zeros(length(Nr), length(Np));


for p = Np
      for r = Nr
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
%                
                    grad_Phi_rhop = sqrt(Nup) .* cos(pm .* phir_) .* besselj_der(pm, beta_rhop .* rhor_) .* beta_rhop;
                    grad_Phi_phip = (-1./rhor_) .* sqrt(Nup) .* pm .* sin(pm .* phir_) .* besselj(pm,  beta_rhop .* rhor_);
                
%                
                
%                
                    
                    grad_Phi_rhor = sqrt(Nur) .* cos(rm .* phir_) .* besselj_der(rm, beta_rhor.* rhor_) .* beta_rhor;
                    grad_Phi_phir = (-1./rhor_) .* sqrt(Nur) .* rm .* sin(rm .* phir_) .* besselj(rm,  beta_rhor .* rhor_);
%                
                
    if (modep  == "TE" && moder  == "TE") || (modep  == "TM" && moder  == "TM")
                    
                    if pm == rm
                            A = Lommel(0, rr, beta_rhop, beta_rhor, pm - 1, rm - 1);
                      
                            D = Lommel(0, rr, beta_rhop, beta_rhor, pm + 1, rm + 1);

                    
                            Icos = intphicos(0, 2*pi, pm, rm);
                            Isin = intphisin(0, 2*pi, pm, rm);
                            K = beta_rhop .* beta_rhor./4;
                      
                            X_til_pr =  sqrt(Nup) .* sqrt(Nur) ...
                            .* K .* (A + D) .* (Icos + Isin);

                            X_til = X_til_pr; 
                      
                    else 
                        
                            X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
                        .* rhor_ .* drho .* dphi;
                            X_til = sum(sum(X_til_pr));
                    end
                    
     elseif (modep == "TE" && moder == "TM")
         
                    X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_phip)...
                       .* rhor_ .* drho .* dphi;
                   X_til = sum(sum(X_til_pr));
                    
     elseif (modep == "TM" && moder == "TE")
         
                   X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_phip)...
                       .* rhor_ .* drho .* dphi;
                   X_til = sum(sum(X_til_pr));
                        
     end             
                    
                    
                    
        
      end
      
end
    