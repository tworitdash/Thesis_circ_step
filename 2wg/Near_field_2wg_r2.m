clear;

F = 20e9;

err = 1; erp = 1; murr = 1; murp = 1;

er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability

epsilonp = erp * er0;   % Permittivity in the medium
mup = mu0 * murp;       % Permeability in the medium
epsilonr = err * er0;   % Permittivity in the medium
mur = mu0 * murr;

rr = 0.0405319403216/1.9;
rp = 0.0405319403216/1.9 * 2;


[fc_1] = fc(rr, err, murr);
[fc_2] = fc(rp, erp, murp);

Nr = find(fc_1 < F);
Np = find(fc_2 < F);

[X_] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr); 

[Spp, Spr, Srp, Srr] = GSM(Nr, Np, F, rp, rr, erp, murp, err, murr, X_);

[rho, phi] = meshgrid(eps:rp/100:rp, -pi-eps:pi/180:pi+eps);

z = -0.02;

[Er] = E_r(Nr, rho, phi, F, rr, z, epsilonr, mur);

Gamma_sum = sum(abs(Srp).^2, 2);

E_aperture = zeros(size(rho));

for k = 1:length(Nr)
    E_aperture = E_aperture + squeeze(Er(k, :, :)) .* Gamma_sum(k);
end

x = rho .* cos(phi);
y = rho .* sin(phi);

figure;
surface(x, y, db(abs(E_aperture))); shading flat;




