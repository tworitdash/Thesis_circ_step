clear;

% M = 10; % Number of elements in between the first and last waveguide

% r = linspace(0.02, 0.03, 10);

L = 0.02; % Length of each waveguide section

F = 4e9:0.5e9:21e9; % Frequency of operation

rt = 0.0405319403216/1.9;
rp = 0.0405319403216/2; % radius of the waveguide
rr = 0.0405319403216/2.1;
rd = 2.2e-2;

r = [rr rp rt];

%% Inner Cross Product of First junction (smaller dimension)

% Nr = 1:1:5; % number of modes on R waveguide
% Np = 1:1:5; % number of modes on P waveguide

N1 = 1:1:5; % number of modes on 1st waveguide
N2 = 1:1:5; % number of modes on 2nd waveguide

erp = 1;
err = 1;
murp = 1;
murr = 1;

% [X_til_rp] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr);
[X_til_rp] = Inner_p(N1, N2, r(2), r(1), erp, murp, err, murr);

%% Inner Cross Product of Last junction (larger dimension)


% Np = 1:1:30; % number of modes on R waveguide
% Nt = 1:1:30; % number of modes on P waveguide

Ns = 1:1:5; % number of modes on last but one waveguide
Ne = 1:1:5; % number of modes on last waveguide

ert = 1;
erp = 1;
murt = 1;
murp = 1;

% [X_til_pt] = Inner_p(Np, Nt, rt, rp, ert, murt, erp, murp);

[X_til_pt] = Inner_p(Ns, Ne, r(end), r(end - 1), ert, murt, erp, murp);

%% Frequency dependent terms:
for k = 1:length(F)

disp('Iteration: ');
disp(k);

% [S33, S34, S43, S44] = GSM(Nr, Np, F(k), rp, rr, erp, murp, err, murr, X_til_rp);
% [S11, S12, S21, S22] = GSM(Np, Nt, F(k), rt, rp, erp, murp, err, murr, X_til_pt);


[S33, S34, S43, S44] = GSM(N1, N2, F(k), r(2), r(1), erp, murp, err, murr, X_til_rp);
[S11, S12, S21, S22] = GSM(Ns, Ne, F(k), r(end), r(end - 1), erp, murp, err, murr, X_til_pt);


Sl = SL(rp, F(k), Ns, L);

   

%   Sl = Sl * SL(rp, F(k), Np, L);
%  [Slr] = SL(rr, F(k), Nr, L);
%  [Slt] = SL(rt, F(k), Nt, L);
 
 [Slr] = SL(r(1), F(k), N1, 0.001);
 [Slt] = SL(r(end), F(k), Ne, 0.001);
 
% Sl = Slr * Slp * Slt;

% I = eye(length(Np), length(Np));

I = eye(length(N2), length(N2));

U1 = inv(I - S22 * Sl * S33 * Sl);
U2 = inv(I - S33 * Sl * S22 * Sl);

STT(k, :, :) = Slt * (S11 + S12 * Sl * U2 * S33 * Sl * S21) * Slt;
STR(k, :, :) = Slt * (S12 * Sl * U2 * S34) * Slr;
SRT(k, :, :) = Slr * (S43 * Sl * U1 * S21) * Slt;
SRR(k, :, :) = Slr * (S44 + S43 * Sl * U1 * S22 * Sl * S34) * Slr;


end
%% Plots


save('Stt3_ratio_1_modes_5_2cm', 'STT');
save('Str3_ratio_1_modes_5_2cm', 'STR');
save('Srt3_ratio_1_modes_5_2cm', 'SRT');
save('Srr3_ratio_1_modes_5_2cm', 'SRR');

