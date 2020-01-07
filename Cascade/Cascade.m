clear;

M = 1; % Number of elements in between the first and last waveguide

L = 0.001; % Length of each waveguide section

F = 4e9:0.5e9:50e9; % Frequency of operation

rt = 0.0405319403216/1.9;
rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;

%% Inner Cross Product of First junction (smaller dimension)

Nr = 1:1:5; % number of modes on R waveguide
Np = 1:1:5; % number of modes on P waveguide

erp = 1;
err = 1;
murp = 1;
murr = 1;

[X_til_rp] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr);

%% Inner Cross Product of Last junction (larger dimension)


Np = 1:1:5; % number of modes on R waveguide
Nt = 1:1:5; % number of modes on P waveguide

ert = 1;
erp = 1;
murt = 1;
murp = 1;

[X_til_pt] = Inner_p(Np, Nt, rt, rp, ert, murt, erp, murp);

%% Frequency dependent terms:
for k = 1:length(F)

disp('Iteration: ');
disp(k);

[S33, S34, S43, S44] = GSM(Nr, Np, F(k), rp, rr, erp, murp, err, murr, X_til_rp);
[S11, S12, S21, S22] = GSM(Np, Nt, F(k), rt, rp, erp, murp, err, murr, X_til_pt);

[Sl] = SL(rp, F(k), Np, L);
[Slr] = SL(rr, F(k), Nr, 0.001);
[Slt] = SL(rt, F(k), Nt, 0.001);

% Sl = Slr * Slp * Slt;

I = eye(length(Np), length(Np));

U1 = inv(I - S22 * Sl * S33 * Sl);
U2 = inv(I - S33 * Sl * S22 * Sl);

STT(k, :, :) = Slt * (S11 + S12 * Sl * U2 * S33 * Sl * S21) * Slt;
STR(k, :, :) = Slt * (S12 * Sl * U2 * S34) * Slr;
SRT(k, :, :) = Slr * (S43 * Sl * U1 * S21) * Slt;
SRR(k, :, :) = Slr * (S44 + S43 * Sl * U2 * S22 * Sl * S34) * Slr;


end
%% Plots


save('Stt3_ratio_1_modes_5_V3', 'STT');
save('Str3_ratio_1_modes_5_V3', 'STR');
save('Srt3_ratio_1_modes_5_V3', 'SRT');
save('Srr3_ratio_1_modes_5_V3', 'SRR');


