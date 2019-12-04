%% Normalization matrix with different modes on the waveguide:
clear;
close all;
c0 = 3e8;    
m = 1; % first digit of the mode number
N = 1:1:50; % second digit of the mode number
%N = 3;
mode = "TM"; % Waveguide mode polarization
% mode = "TM"

%<<<<<<< HEAD
% F =    1.4132e+11;
%=======
% F = 1e9:1e9:20e9;
F = 200e9;
%>>>>>>> 054067d7bb6b39fdb9b3e302b4d131398d191f0f

r = 0.0405319403216/2; % radius of the waveguide
er = 1; % relative  permittivity
mur = 1; % relative Permeability
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;

drho = r/1000;
dphi = pi/2000;

[rho_, phi_] = meshgrid(eps:drho:r, eps:dphi:2*pi-eps);  % domain for the fields on one cross-section of the waveguide

z = 0; 
% Q = zeros(size(F, 2), N(end), N(end));
Q = zeros(N(end), N(end));
Q_numerical = zeros(N(end), N(end));
Z = zeros(1, size(F, 2));
Y = zeros(1, size(F, 2));

for k = 1:length(F)

for i = 1:length(N)
% <<<<<<< HEAD
    disp(N(i));
    
    omega = 2 * pi * F(k);

    lamb = c0./F(k); % wavelength
 
    beta = 2 * pi ./ lamb;
% =======
% >>>>>>> 054067d7bb6b39fdb9b3e302b4d131398d191f0f
    [Erho, Ephi, Ez, Hrho, Hphi, Hz, beta_z] = E_and_H(rho_, phi_, er, mur, z, r, m, N(i), mode, F(k));
    beta_rho = sqrt((beta.^2 - beta_z.^2));
%% Analytical    
    if mode == "TE"
        K = beta_z .* beta_rho.^2 ./ (4 .* omega .* mu * epsilon.^2);
    elseif mode == "TM"
        K = beta_z .* beta_rho.^2 ./ (4 .* omega .* mu.^2 * epsilon);
    end
    
    A = Lommel(0, r, beta_rho, beta_rho, m - 1, m - 1);
    B = Lommel(0, r, beta_rho, beta_rho, m - 1, m + 1);
    C = Lommel(0, r, beta_rho, beta_rho, m + 1, m + 1);
    Isin = intphisin(0, 2*pi, m, m);
    Icos = intphicos(0, 2*pi, m, m);
    
    Qij = K .* ((Isin + Icos) .* (A + C) + 2 .* B .* (Isin - Icos));
    
%     Poyn = (Erho .* Hphi - Hrho .* Ephi) .* rho_ * drho .* dphi;
%     Qij = sum(sum(Poyn));
%     
    Q(i, i) = Qij;
%% Numerical Q
    Poyn = (Erho .* Hphi - Hrho .* Ephi) .* rho_ * drho .* dphi;
    Qij_numerical = sum(sum(Poyn));
    
    Q_numerical(i, i) = Qij_numerical;
%     if mode == "TE"
%         Z_i = 2 * pi * F(k) * mu./ beta_z;
%     elseif mode == "TM"
%         Z_i = beta_z ./ (2 * pi * F(k) .* epsilon);
%     end
%     Z(k) = Z_i;
%     Y(k) = 1./Z(k);
end
end

% <<<<<<< HEAD
% figure;
% plot(m, db(abs(diag(Q)))/10, 'LineWidth', 2); grid on;
% 
% xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Normalization Constant Q_{1, n}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Normalization Constant for', mode,'_{1, n} modes'], 'FontSize', 12, 'FontWeight', 'bold');

% figure;
% plot(N, (real(diag(Q))), 'LineWidth', 2); grid on;
% hold on;
% plot(N, (imag(diag(Q))), 'LineWidth', 2); grid on;
% 
% xlabel('n in TE_{m, 2} modes', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Normalization Constant Q_{m, 2}', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Normalization Constant for', mode,'_{m, 2} modes'], 'FontSize', 12, 'FontWeight', 'bold');
% =======
% figure;
% plot(N, db(abs(diag(Q)))/10, 'LineWidth', 2); grid on;
% 
% xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Normalization Constant Q_{1, n}(dB)', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Normalization Constant for', mode,'_{1, n} modes'], 'FontSize', 12, 'FontWeight', 'bold');
% % 
% figure;
% plot(N, real(diag(Q)), 'LineWidth', 2); grid on;
% hold on;
% plot(N, imag(diag(Q)), 'LineWidth', 2); grid on;
% hold on;
% plot(N, real(diag(Q_numerical)), '-.', 'LineWidth', 2); grid on;
% hold on;
% plot(N, imag(diag(Q_numerical)), '-.', 'LineWidth', 2); grid on;
% 
% xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Normalization Constant Q_{1, n}', 'FontSize', 12, 'FontWeight', 'bold');
% title(['Normalization Constant for', mode,'_{1, n} modes'], 'FontSize', 12, 'FontWeight', 'bold');
% % >>>>>>> 054067d7bb6b39fdb9b3e302b4d131398d191f0f
% legend({'Analytical Re(Q)', 'Analytical Im(Q)', 'Numerical Re(Q)', 'Numerical Im(Q)'}, ...
%     'FontSize', 12, 'FontWeight', 'bold');

figure;
plot(N, real(diag(Q)), 'LineWidth', 2); grid on;
hold on;
plot(N, imag(diag(Q)), 'LineWidth', 2); grid on;
hold on;
plot(N, real(diag(Q_numerical)), '*', 'LineWidth', 2); grid on;
hold on;
plot(N, imag(diag(Q_numerical)), '*', 'LineWidth', 2); grid on;

xlabel('n in TE_{1, n} modes', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Normalization Constant Q_{1, n}', 'FontSize', 12, 'FontWeight', 'bold');
title(['Normalization Constant for', mode,'_{1, n} modes'], 'FontSize', 12, 'FontWeight', 'bold');
% >>>>>>> 054067d7bb6b39fdb9b3e302b4d131398d191f0f
legend({'Analytical Re(Q)', 'Analytical Im(Q)', 'Numerical Re(Q)', 'Numerical Im(Q)'}, ...
    'FontSize', 12, 'FontWeight', 'bold');
% hold on;
% hold on;
% plot(N, real(diag(Q))/2, 'LineWidth', 2); grid on;
% hold on;
% plot(N, imag(diag(Q)/2), 'LineWidth', 2); grid on;
% hold on;

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

% figure;
% plot(F * 1e-9, real(Y), 'LineWidth', 2); grid on;
% hold on;
% plot(F * 1e-9, imag(Y), 'LineWidth', 2);
% 
% xlabel('Frequency (GHz)', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('Admittance Y (Mho)', 'FontSize', 12, 'FontWeight', 'bold');
% title('Wave Admittance', 'FontSize', 12, 'FontWeight', 'bold');
% legend({'Re(Y)', 'Im(Y)'}, 'FontSize', 12, 'FontWeight', 'bold');
