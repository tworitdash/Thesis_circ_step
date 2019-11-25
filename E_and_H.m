
<<<<<<< HEAD
function [Erho, Ephi, Ez, Hrho, Hphi, Hz, beta_z, xmn] = E_and_H(rho_, phi_, er, mur, z, r, m, n, mode, f)
=======
function [Erho, Ephi, Ez, Hrho, Hphi, Hz, beta_z] = E_and_H(rho_, phi_, er, mur, z, r, m, n, mode, f)
>>>>>>> 054067d7bb6b39fdb9b3e302b4d131398d191f0f


%% Waveguide parameters 


c0 = 3e8;             % Speed of light in vacuum
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;       % Permeability in the medium
%% Zeros of Bessel's functions and their derivatives:

xm = linspace(0.1, 10000, 100000);

Jm = @(z) besselj(m, z);                                         % Bessel's function

dz = 1e-5;

Jm_der = @(z) (besselj(m, z + dz) - besselj(m, z - dz))./(2 * dz); % Derivative of the Bessel's function

if mode == "TE"
    fun = Jm_der;
else
    fun = Jm;
end

ym = fun(xm);                   % Values of the Bessel's function at the positions defined by xm

chsign = find(diff(sign(ym)));  % Detection of the sign changes 

xmn_all = zeros(size(chsign));

for i = 1:size(chsign, 2)
    xmn_all(i) = fzero(fun, xm(chsign(i)));  % finding the roots near the points where Jm changes sign
end



xmn = xmn_all(n); 


beta_rho = xmn/r;  % wave number along the rho direction (\beta_{\rho})

%% Cut of frequency

fc = xmn ./ (2 * pi * r * sqrt(mu .* epsilon));
<<<<<<< HEAD
disp(fc);
=======
>>>>>>> 054067d7bb6b39fdb9b3e302b4d131398d191f0f

%% Source Parameters
% f = fc + eps;

omega = 2 * pi * f;

lamb = c0./f; % wavelength
 
beta = 2 * pi ./ lamb; % wave number

beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));

%[rho_, phi_] = meshgrid(rho, phi);


%% Wave equations

if mode == "TE"
    
    [Erho, Ephi, Ez] = E_TE(epsilon, m, rho_, phi_, beta_rho, z, beta);
    [Hrho, Hphi, Hz] = H_TE(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
   
elseif mode == "TM"
    [Erho, Ephi, Ez] = E_TM(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
    [Hrho, Hphi, Hz] = H_TM(mu, m, rho_, phi_, beta_rho, z, beta);
end


%% Cartesian Coordinates 

% Ex = cos(phi_) .* Erho - sin(phi_) .* Ephi;
% Ey = sin(phi_) .* Erho + cos(phi_) .* Ephi;
% 
% x = rho_ .* cos(phi_);
% y = rho_ .* sin(phi_);


% Hx = cos(phi_) .* Hrho - sin(phi_) .* Hphi;
% Hy = sin(phi_) .* Hrho + cos(phi_) .* Hphi;

% E = sqrt(abs(Erho).^2 + abs(Ephi).^2 + Ez.^2);
% H = sqrt(abs(Hrho).^2 + abs(Hphi).^2 + Hz.^2);


% dx = x(2, 1) - x(1, 1);
% dy = y(1, 2) - y(1, 1);
% 
% 
% dHx = gradient(imag(Hx), dx);
% dHy = gradient(imag(Hy), dy);
% 
% dEx = gradient(real(Ex), dx);
% dEy = gradient(real(Ey), dy);


end
% figure;
% contour(x, y, abs(E));
% 
% grid on;
% hold on;
% 
% quiver(x(1:e:end, 1:e:end), y(1:e:end, 1:e:end), dEx(1:e:end, 1:e:end), dEy(1:e:end, 1:e:end));
% title('E', 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
% contour(x, y, abs(H));
% 
% grid on;
% hold on;
% 
% quiver(x(1:e:end, 1:e:end), y(1:e:end, 1:e:end), dHx(1:e:end, 1:e:end), dHy(1:e:end, 1:e:end));
% title('H', 'FontSize', 12, 'FontWeight', 'bold');

% e = 50;
% 
% figure;
% h = pcolor(x,y, abs(E));
% 
% set(h,'ZData',-1+zeros(size(E)))
% hold on;
% quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),real(Hx(1:e:end, 1:e:end)),real(Hy(1:e:end, 1:e:end)), 2, 'w');
% shading interp;
% colorbar;
% colormap('jet');
% title('E', 'FontSize', 12, 'FontWeight', 'bold');
% 
% 
% figure;
% h = pcolor(x,y, abs(H));
% 
% set(h,'ZData',-1+zeros(size(H)))
% hold on;
% 
% shading interp;
% 
% quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),imag(Ex(1:e:end, 1:e:end)),imag(Ey(1:e:end, 1:e:end)), 2, 'w');
% colorbar;
% colormap('jet');
% title('H', 'FontSize', 12, 'FontWeight', 'bold');

