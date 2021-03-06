clear
%% Parameters of the waveguides

rp = 0.0405319403216/2;   % radius of the bigger waveguide
rr = 0.0405319403216/2.1; % radius of the smaller waveguide

% rp = 0.0405319403216/1.9 * 2;   % radius of the bigger waveguide
% rr = 0.0405319403216/1.9; % radius of the smaller waveguide

erp = 1;                  % Relative Permittivity of P waveguide
err = 1;                  % Relative Permittivity of R waveguide
murp = 1;                 % Relative Permeability of P waveguide
murr = 1;                 % Relative Permeability of R waveguide

r = [rr rp];

F = 4e9:0.5e9:21e9; % Frequency of operation

%% Inner cross product

Nr = 1:1:20; % number of modes on 1st waveguide
Np = 1:1:20; % number of modes on 2nd waveguide



[X_til_rp] = Inner_p(Nr, Np, r(2), r(1), erp, murp, err, murr);

parfor k = 1:length(F)
    
    disp('Frequency Iteration')
    disp(k);
    
    [Spp_, Spr_, Srp_, Srr_] = GSM(Nr, Np, F(k), r(2), r(1), erp, murp, err, murr, X_til_rp); % GSM 
    
    slr = SL(rr, F(k), Nr, 0.001);  % Phase due to the length (height) of the R cylinder
    slp = SL(rp, F(k), Np, 0.001);  % Phase due to the length (height) of the P cylinder
    
%     Spp(k, :, :) = slp * Spp_ * slp';
%     Spr(k, :, :) = slp * Spr_ * slr;
%     Srp(k, :, :) = slr * Srp_ * slp;
%     Srr(k, :, :) = slr * Srr_ * slr';
%     
    
    Spp(k, :, :) =  slp * Spp_ * slp;
    Spr(k, :, :) =  slp * Spr_ * slr;
    Srp(k, :, :) =  slr * Srp_ * slp;
    Srr(k, :, :) =  slr * Srr_ * slr;
%     Srr(k, :, :) = Srr_;
end

%% Saving the Data
save('Stt2_ratio_1_modes_20_fc_align', 'Spp');
save('Str2_ratio_1_modes_20_fc_align', 'Spr');
save('Srt2_ratio_1_modes_20_fc_align', 'Srp');
save('Srr2_ratio_1_modes_20_fc_align', 'Srr');
