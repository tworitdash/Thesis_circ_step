clear;

F = 14e9;

err = 1; erp = 1; murr = 1; murp = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;

% <<<<<<< HEAD
rr = 0.0405319403216/1.9;
rp = 0.0405319403216/1.9 * 2;
% =======
% rr = 0.0405319403216/2.1;
% rp = 0.0405319403216/2;
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608


[fc_1] = fc(rr, err, murr);
[fc_2] = fc(rp, erp, murp);

Nr1 = find(fc_1 < F);
Np1 = find(fc_2 < F);

% <<<<<<< HEAD
Nr = 1:1:10;
Np = 1:1:40;
% =======
% Nr = 1:1:20;
% Np = 1:1:20;
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608
% Np = Nr;

[X_] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr); 

[Spp, Spr, Srp, Srr] = GSM(Nr, Np, F, rp, rr, erp, murp, err, murr, X_);

%% 
[rho, phi] = meshgrid(eps:rr/100:rr, -pi-eps:pi/180:pi+eps);

z = 0.001;

[Er_rho, Er_phi, Er_z] = E_r(Nr1, rho, phi, F, rr, z, epsilonr, mur);

% <<<<<<< HEAD
Gamma_sum = 1 + sum(Srr, 2);
% =======
% Gamma_sum = 1 + sum(Srr.^2, 2);
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:length(Nr1)
    E_aperture_rho = E_aperture_rho + squeeze(Er_rho(k, :, :)) .* Gamma_sum(k);
    E_aperture_phi = E_aperture_phi + squeeze(Er_phi(k, :, :)) .* Gamma_sum(k);
    E_aperture_z = E_aperture_z + squeeze(Er_z(k, :, :)) .* Gamma_sum(k);
end

E_aperture = sqrt(E_aperture_rho.^2 + E_aperture_phi.^2 + E_aperture_z.^2);
% E_aperture = sqrt(E_aperture_rho.^2 + E_aperture_phi.^2);

x = rho .* cos(phi);
y = rho .* sin(phi);
% 
figure;
surface(x, y, ((abs(E_aperture))./max(abs(E_aperture)))); shading flat;
% colormap('jet');



[ficname,pathname] = uigetfile('*.efe','fichier ''.efe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');
% <<<<<<< HEAD
% system(['"C:\Program Files (x86)\GnuWin32\bin\sed" -e "/^#/d;/^*/d" ',' "',nomfic,'"| "C:\Program Files (x86)\GnuWin32\bin\tr" -s " " " " > result.txt']);  
% =======
system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608
A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

% <<<<<<< HEAD
% [rho_f, phi_f] = meshgrid(rho_f, pi - phi_f);
% =======
[rho_f, phi_f] = meshgrid(rho_f, phi_f);
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608

x_f = rho_f .* cos(phi_f);
y_f = rho_f .* sin(phi_f);

f = 20*36000:1:21*36000 - 1;

E_rho = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_phi = sqrt(aux(f, 3).^2 + aux(f, 4).^2);
E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2 + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);

figure;

surface(x_f, y_f, ((abs(E_tot_reshape'))./max(abs(E_tot_reshape')))); shading flat;


% <<<<<<< HEAD
figure;

plot(phi(:, 1), db(abs(E_aperture(:, 30)/max(abs(E_aperture(:, 30))))));
hold on;

plot(phi_f(:, 1), db(abs(E_tot_reshape(30, :)'/max(abs(E_tot_reshape(30, :)))')));

grid on;


%----------------------------------------------------------------------------------------------------------


% =======
% %%
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608
[rho, phi] = meshgrid(eps:rp/100:rp, -pi-eps:pi/180:pi+eps);

z = -0.001;

[Ep_rho, Ep_phi, Ep_z] = E_r(Np1, rho, phi, F, rp, z, epsilonp, mup);

% <<<<<<< HEAD
Gamma_sum = sum(Spr, 2);% + sum(Spp, 2);
% =======
% Gamma_sum = sum(Spr.^2, 2);
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608

E_aperture_rho = zeros(size(rho));
E_aperture_phi = zeros(size(rho));
E_aperture_z = zeros(size(rho));

for k = 1:length(Np1)
    E_aperture_rho = E_aperture_rho + squeeze(Ep_rho(k, :, :)) .* (Gamma_sum(k));
    E_aperture_phi = E_aperture_phi + squeeze(Ep_phi(k, :, :)) .* (Gamma_sum(k));
    E_aperture_z = E_aperture_z + squeeze(Ep_z(k, :, :)) .* (Gamma_sum(k));
end

E_aperture = sqrt(E_aperture_rho.^2 + E_aperture_phi.^2 + E_aperture_z.^2);

% E_aperture = sqrt(E_aperture_rho.^2 + E_aperture_phi.^2);

x = rho .* cos(phi);
y = rho .* sin(phi);
% 
figure;
surface(x, y, db((abs(E_aperture))./max(abs(E_aperture)))); shading flat;



[ficname,pathname] = uigetfile('*.efe','fichier ''.efe'' a convertir ?');
nomfic = [pathname ficname];
i0 = find(ficname=='.');
% <<<<<<< HEAD
% system(['"C:\Program Files (x86)\GnuWin32\bin\sed" -e "/^#/d;/^*/d" ',' "',nomfic,'"| "C:\Program Files (x86)\GnuWin32\bin\tr" -s " " " " > result.txt']);    
% =======
system(['sed -e "/^#/d;/^*/d" ',' "',nomfic,'"| tr -s " " " " > result.txt']);  
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608
A = load('result.txt');

rho_f = A(1:100, 1);
phi_f = A(1:100:36000, 2) * pi/180;
z_f = A(1, 1);
aux = A(:, 4:9);

% <<<<<<< HEAD
[rho_f, phi_f] = meshgrid(rho_f, pi - phi_f);
% =======
% [rho_f, phi_f] = meshgrid(rho_f, phi_f);
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608

x_f = rho_f .* cos(phi_f);
y_f = rho_f .* sin(phi_f);

f = 20*36000:1:21*36000 - 1;

E_rho = sqrt(aux(f, 1).^2 + aux(f, 2).^2);
E_phi = sqrt(aux(f, 3).^2 + aux(f, 4).^2);
E_z = sqrt(aux(f, 5).^2 + aux(f, 6).^2);

E_tot = sqrt(abs(E_rho).^2 + abs(E_phi).^2 + abs(E_z).^2);

E_tot_reshape = reshape(E_tot, 100, 360);

figure;

surface(x_f, y_f, db((abs(E_tot_reshape'))./max(abs(E_tot_reshape')))); shading flat;

% <<<<<<< HEAD
figure;

plot(phi(:, 1), db(abs(E_aperture(:, 30)/max(abs(E_aperture(:, 30))))));
hold on;

plot(phi_f(:, 1), db(abs(E_tot_reshape(30, :)'/max(abs(E_tot_reshape(30, :)))')));

grid on;

% =======
% %% 
% 
% 
% 
% % colormap('jet');
% figure;
% 
% plot(rho(1, :), db(abs(E_aperture(30, :)/max(abs(E_aperture(30, :))))));
% hold on;
% 
% plot(rho_f(1, :), db(abs(E_tot_reshape(:, 40)'/max(abs(E_tot_reshape(:, 40)))')));
% >>>>>>> 37893f8e2e343fb0a1d771a64115f47cd4bdd608

