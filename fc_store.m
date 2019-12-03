
r = 0.0405319403216/2;

er = 1;
mur = 1;
c0 = 3e8;             % Speed of light in vacuum
er0 = 8.85418782e-12; % Free space permittivity
mu0 = 1.25663706e-6;  % Free Space Permeability
epsilon = er * er0;   % Permittivity in the medium
mu = mu0 * mur; 

x_TE = csvread('TE_Xmn');
x_TM = csvread('TM_Xmn');

fc_TE = x_TE ./ (2 * pi * r * sqrt(mu .* epsilon));

fc_TM = x_TM ./ (2 * pi * r * sqrt(mu .* epsilon));

save('fc_TE', 'fc_TE');
save('fc_TM', 'fc_TM');

%% 
c_ = load('fc_TE');
d_ = c_.fc_TE;
d_(1, 2)