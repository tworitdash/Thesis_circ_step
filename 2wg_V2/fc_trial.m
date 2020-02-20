er = 1;
mur = 1;


r = 0.0405319403216/2.1;

c0 = 3e8;
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur;

Str = load('Xmn.mat');
Xmn = Str.Xmn;

for i = 1:size(Xmn, 2)
    X(i) = Xmn(i).xmn;
    M(1, i) = Xmn(i).mode;
    M(2, i) = Xmn(i).m;
    M(3, i) = Xmn(i).n;
end


fc_ = X ./ (2 * pi * r * sqrt(mu .* epsilon));