%% Normalization matrix with different modes on the waveguide:
function [Q, Z, Y, K] = QZcalculation_v2(N, F, r, er, mur, rho_, phi_, z, drho, dphi)
c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;


% Q = zeros(size(F, 2), N(end), N(end));
Q = zeros(length(N), length(N));
Z = zeros(length(N), length(N));
Y = zeros(length(N), length(N));


xmn_ = zeros(length(N));
% beta_z_ = zeros(N);

K = zeros(length(N), length(N));


%Y = zeros(m(end), m(end));

omega = 2 * pi * F;

lamb = c0./F; % wavelength
 
beta = 2 * pi ./ lamb;



for i = 1:length(N)
    
    disp(i);
 %% Analytical Q
 
       Str = load('Xmn.mat');
       Xmn = Str.Xmn;
      
      
      xmn_(i) = Xmn(N(i)).xmn;
      m = Xmn(N(i)).m;
      
      mode = Xmn(N(i)).mode;
      
      fc = xmn_(i) ./ (2 * pi * r * sqrt(mu .* epsilon));
      disp(fc);
      
      beta_rho = xmn_(i)./r;

      beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
%         beta_z = sqrt((beta.^2 - beta_rho.^2));
    
      if mode == "TE"

        K(i, i) = beta_z ./ (omega .* mu .* epsilon^2);
        
      elseif mode == "TM"

        K(i, i) = beta_z ./ (omega .* mu.^2 .* epsilon);
        
      end
%       
%                 if mode == "TE"
%                     Norm = (epsilon * pi/2 .* (xmn_(i).^2 - m.^2) .* (besselj(m, xmn_(i))).^2).^(-1);
%                 elseif mode == "TM"
%                     Norm = (epsilon .* pi/2 .* xmn_(i).^2 .* (besselj_der(m, xmn_(i))).^2).^(-1);
%                 end
        Norm = 1;        
               
      Const = K(i, i) .* beta_rho.^2 ./ 4 .* Norm;  

        A = Lommel(0, r, beta_rho, beta_rho, m - 1, m - 1);

        C = Lommel(0, r, beta_rho, beta_rho, m + 1, m + 1);
        
        if m == 0
            Isin = 0;
            Icos = 2 * pi;
        else
            Isin = intphisin(0, 2*pi, m, m);
            Icos = intphicos(0, 2*pi, m, m);
        end
    
       Qij = Const .* ((Isin + Icos) .* (A + C));
    
       Q(i, i) = Qij;
      
%      
%% Impedance (Z) and Admittance (Y) 
      
    
    if mode == "TE"
        Z_i = 2 * pi * F * mu./ beta_z;
        Y_i = 1./Z_i;
    elseif mode == "TM"
        Z_i = beta_z ./ (2 * pi * F .* epsilon);
        Y_i = 1./Z_i;
    end
    
%     Z(i , i) = Z_i;
%     Y(i, i) = Y_i;
    
    Z(i, i) = Z_i;
    Y(i, i) = Y_i;
    
    
%     beta_rho = xmn_./r;

%     beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
  
    
end


end

% figure;
% plot(m, db(abs(diag(Q)))/10, 'LineWidth', 2); grid on;
% 
% xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Normalization Constant Q_{1, n}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Normalization Constant for', mode,'_{1, n} modes'], 'FontSize', 12, 'FontWeight', 'bold');

% figure;
% plot(m, (real(diag(Q))), 'LineWidth', 2); grid on;
% hold on;
% plot(m, (imag(diag(Q))), 'LineWidth', 2); grid on;
% 
% xlabel('n in TE_{m, 2} modes', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Normalization Constant Q_{m, 2}', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Normalization Constant for', mode,'_{m, 2} modes'], 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Q)', 'Im(Q)'}, 'FontSize', 12, 'FontWeight', 'bold');

% figure;
% 
% plot(F * 1e-9, real(Z), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, imag(Z), 'LineWidth', 2);
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Impedance Z (\Omega)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Wave Impedance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Z)', 'Im(Z)'}, 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
% plot(F * 1e-9, real(Y), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, imag(Y), 'LineWidth', 2);
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Admittance Y (Mho)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Wave Admittance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Y)', 'Im(Y)'}, 'FontSize', 12, 'FontWeight', 'bold');