close all;
clear;
clc;

%% Waveguide parameters
d = 0.0405319403216; % diameter in meters
r = d/2;             % radius


c0 = 3e8;
er0 = 8.85418782e-12;
mu0 = 1.25663706e-6;

er = 1;                 %relative permittivity
mur = 1;                %relative permeability
epsilon = er * er0;
mu = mu0 * mur;
%% Zeros of Bessel's functions and their derivatives:
m = input('Input the first digit of the mode number: ');
n = input('Input the second digit of the mode number: ');
mode = input('Which mode (TE/TM): ', 's');

xm = linspace(0.1, 100, 100000);

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

z = 0.2;

%% Cut of frequency

fc = xmn ./ (2 * pi * r * sqrt(mu .* epsilon));

%% Source Parameters
f = fc + eps;

omega = 2 * pi * f;

lamb = c0./f; % wavelength
 
beta = 2 * pi ./ lamb; % wave number

rho = eps:r/1000:r;

phi = eps:pi/2000:2*pi-eps;

[rho_, phi_] = meshgrid(rho, phi);


%% Wave equations

if mode == "TE"
    
    [Erho, Ephi, Ez] = E_TE(epsilon, m, rho_, phi_, beta_rho, z, beta);
    [Hrho, Hphi, Hz] = H_TE(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
   
elseif mode == "TM"
    [Erho, Ephi, Ez] = E_TM(epsilon, m, rho_, phi_, beta_rho, z, beta, omega, mu);
    [Hrho, Hphi, Hz] = H_TM(mu, m, rho_, phi_, beta_rho, z, beta);
end


%% Cartesian Coordinates 

Ex = cos(phi_) .* Erho - sin(phi_) .* Ephi;
Ey = sin(phi_) .* Erho + cos(phi_) .* Ephi;

x = rho_ .* cos(phi_);
y = rho_ .* sin(phi_);


Hx = cos(phi_) .* Hrho - sin(phi_) .* Hphi;
Hy = sin(phi_) .* Hrho + cos(phi_) .* Hphi;

% figure;
% quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),abs(Hx(1:e:end, 1:e:end)),abs(Hy(1:e:end, 1:e:end)));
% 
% figure;
% quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),abs(Ex(1:e:end, 1:e:end)),abs(Ey(1:e:end, 1:e:end)));

% figure;
% surface(x, y, abs(sqrt(abs(Ex).^2 + abs(Ey).^2))); shading flat;
% %colormap('jet');
% title('E', 'FontSize', 12, 'FontWeight', 'bold');
% colorbar;
% hold on;
% quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),abs(Ex(1:e:end, 1:e:end)),abs(Ey(1:e:end, 1:e:end)), 0.7, 'w');
% 
% figure;
% surface(x, y, abs(sqrt(abs(Hx).^2 + abs(Hy).^2 + abs(Hz).^2))); shading flat;
% %colormap('jet');
% title('H', 'FontSize', 12, 'FontWeight', 'bold');
% colorbar;
% hold on;
% quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),abs(Hx(1:e:end, 1:e:end)),abs(Hy(1:e:end, 1:e:end)), 0.7, 'r');

E = sqrt(abs(Erho).^2 + abs(Ephi).^2 + Ez.^2);
H = sqrt(abs(Hrho).^2 + abs(Hphi).^2 + Hz.^2);


dx = x(2, 1) - x(1, 1);
dy = y(1, 2) - y(1, 1);


dHx = gradient(imag(Hx), dx);
dHy = gradient(imag(Hy), dy);

dEx = gradient(real(Ex), dx);
dEy = gradient(real(Ey), dy);

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

e = 100;

% figure;
% h = pcolor(x,y, abs(E));
% 
% set(h,'ZData',-1+zeros(size(E)))
% hold on;
% quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),imag(Hx(1:e:end, 1:e:end)),imag(Hy(1:e:end, 1:e:end)), 2, 'w');
% shading interp;
% colorbar;
% colormap('jet');
% xlabel('\rho [m]', 'FontSize', 12, 'FontWeight', 'bold');
% ylabel('\rho [m]', 'FontSize', 12, 'FontWeight', 'bold');
% title([mode, '_{', num2str(m), num2str(n), '} Linear Scale'], 'FontSize', 12, 'FontWeight', 'bold');
% print([mode, '_', num2str(m), num2str(n)], '-depsc');


figure;
h = pcolor(x,y, abs(H));

set(h,'ZData',-1+zeros(size(H)))
hold on;

shading interp;

quiver(x(1:e:end, 1:e:end),y(1:e:end, 1:e:end),imag(Ex(1:e:end, 1:e:end)),imag(Ey(1:e:end, 1:e:end)), 2, 'w');
colorbar;
colormap('jet');

xlabel('\rho [m]', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('\rho [m]', 'FontSize', 12, 'FontWeight', 'bold');
title([mode, '_{', num2str(m), num2str(n), '} Linear Scale'], 'FontSize', 12, 'FontWeight', 'bold');
print([mode, '_', num2str(m), num2str(n)], '-depsc');

% 
% figure;
%    
% polarPcolor(rho, phi.*180/pi, abs(Erho)); shading flat;
% colormap('jet');
% 
% title('E_{\rho}', 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
%    
% polarPcolor(rho, phi.*180/pi, abs(Ephi)); shading flat;
% colormap('jet');
% title('E_{\phi}', 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
% polarPcolor(rho, phi.*180/pi, abs(Hrho)); shading flat;
% colormap('jet');
% title('H_{\rho}', 'FontSize', 12, 'FontWeight', 'bold');
% 
% figure;
%    
% polarPcolor(rho, phi.*180/pi, abs(Hphi)); shading flat;
% colormap('jet');
% title('H_{\phi}', 'FontSize', 12, 'FontWeight', 'bold');

% 
% figure;
% polarPcolor(rho, phi.*180/pi, abs(Hz)); shading flat; colormap('jet');
% title('Hz', 'FontSize', 12, 'FontWeight', 'bold');
% colormap('jet');
% figure;
% surface(rho, phi.*180/pi, abs(Erho)); shading flat; colormap('jet');
% 
% figure;
% surface(rho, phi.*180/pi, abs(Hphi)); shading flat; colormap('jet');



% figure;
% surface(rho, phi.*180/pi, abs(E)); shading flat; colormap('jet');
% 
% figure;
% surface(rho, phi.*180/pi, abs(H)); shading flat; colormap('jet');


% figure;
% polarPcolor(rho, phi.*180/pi, abs(E)); shading flat; colormap('jet');
% title('E', 'FontSize', 12, 'FontWeight', 'bold');
% colormap('jet');
% 
% 
% figure;
% polarPcolor(rho, phi.*180/pi, abs(H)); shading flat; colormap('jet');
% title('H', 'FontSize', 12, 'FontWeight', 'bold');
% colormap('jet');

% figure;
% contour(rho, phi.*180/pi, abs(E));
% figure;
% contour(rho, phi.*180/pi, abs(H));
% figure;
% polarPcolor(rho, phi.*180/pi, abs(sqrt(abs(Hrho).^2 + abs(Hphi).^2))); shading flat; colormap('jet');
% title('H_{\rho} + H_{\phi}', 'FontSize', 12, 'FontWeight', 'bold');



% figure;
% contour(rho, phi.* 180/pi, abs(E));
%    
% figure;
%    
% polarPcolor(rho, phi.*180/pi, abs(H)); shading flat;
% 
% figure;
% contour(rho, phi.* 180/pi, abs(H));
