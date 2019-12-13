%% Normalization matrix with different modes on the waveguide:
function [Q, Z, Y, xmn_] = QZcalculation(m, N, mode, F, r, er, mur, rho_, phi_, z, drho, dphi)

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;


% Q = zeros(size(F, 2), N(end), N(end));
Q = zeros(N(end), N(end));
Z = zeros(N(end), N(end));
Y = zeros(N(end), N(end));
xmn_ = zeros(1, N(end));

%Y = zeros(m(end), m(end));



for i = 1:length(N)
    disp(N(i));
    
    [Erho, Ephi, Ez, Hrho, Hphi, Hz, beta_z, xmn_i] = E_and_H(rho_, phi_, er, mur, z, r, m, N(i), mode, F);
    
    xmn_(i, i) = xmn_i;
    
    Poyn = (Erho .* Hphi - Hrho .* Ephi) .* rho_ * drho .* dphi;
    Qij = sum(sum(Poyn));
    
    Q(i, i) = Qij;
    
    if mode == "TE"
        Z_i = 2 * pi * F * mu./ beta_z;
        Y_i = 1./Z_i;
    elseif mode == "TM"
        Z_i = beta_z ./ (2 * pi * F .* epsilon);
        Y_i = 1./Z_i;
    end
    
    Z(i , i) = Z_i;
    Y(i, i) = Y_i;
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
