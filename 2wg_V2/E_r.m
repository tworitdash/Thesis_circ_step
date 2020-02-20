function [Erho, Ephi, Ez] = E_r(Nr, rho, phi, F, r, z, epsilon, mu)
        c0 = 3e8;
        Str = load('Xmn.mat');
        Xmn = Str.Xmn;
        
        E_i = zeros(length(Nr), size(rho, 1), size(rho, 2));
        
        for i = 1:length(Nr)
            disp(i);
            mode = Xmn(i).mode;
            m = Xmn(i).m;
            xmn = Xmn(i).xmn;
            
            beta_rho = xmn./r;
            beta = (2 .* pi .* F)./c0;
            
            omega = 2 .* pi .* F;
            
            if mode == "TE"
                    A = 1;
   
                    C1 = 1; % C1 and D1 are for the rho component
                    D1 = 0; 

                    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));

                    Erho(i, :, :) = -A .* m./(epsilon .* rho) .* besselj(m, beta_rho .* rho) .* (-C1 .* sin(m .* phi)...
                        + D1 .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
                    Ephi(i, :, :) = A .* beta_rho./epsilon .* besselj_der(m, beta_rho .* rho) .* (C1 .* cos(m .* phi)...
                        + D1 .* sin(m .* phi)) .* exp(-1j .* beta_z .* z);
                    Ez(i, :, :) = zeros(size(rho));
                    
%                     E_i(i, :, :) = sqrt(Erho.^2 + Ephi.^2 + Ez.^2);
                
            else
                    B = 1;

                    C = 1;
                    D = 0;

                    beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));

                    Erho(i, :, :)  = -B .* (beta_rho .* beta_z ./ (omega .* mu .* epsilon)) .* besselj_der(m, beta_rho .* rho).* (C .* cos(m .* phi)...
                        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
                    Ephi(i, :, :)  = -B .* (m .* beta_z ./ (omega .* mu .* epsilon .* rho)) .* besselj(m, beta_rho .* rho) .* (- C .* sin(m .* phi)...
                        + D .* cos(m .* phi)) .* exp(-1j .* beta_z .* z);
                    Ez(i, :, :)  = -1j .* B .* (beta_rho.^2 ./ (omega .* mu .* epsilon)) .* besselj(m, beta_rho .* rho) .*  (C .* cos(m .* phi)...
                        + D .* sin(m .* phi)) * exp(-1j .* beta_z .* z);
                    
%                     E_i(i, :, :) = sqrt(Erho.^2 + Ephi.^2 + Ez.^2);
            end
            
            
        end
        
%         E = E_i;
    
end