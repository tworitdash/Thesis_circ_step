function [STT, STR, SRT, SRR] = cascade_3(r, L, F, N, err, murr, erp, murp, ert, murt)

%% Inner cross products of the first and last junction


[X_til_rp] = Inner_p(N, N, r(2), r(1), erp, murp, err, murr);

[X_til_pt] = Inner_p(N, N, r(3), r(2), ert, murt, erp, murp);

%% Frequency dependent terms:
for k = 1:length(F)

disp('Iteration: ');
disp(k);

% [S33, S34, S43, S44] = GSM(Nr, Np, F(k), rp, rr, erp, murp, err, murr, X_til_rp);
% [S11, S12, S21, S22] = GSM(Np, Nt, F(k), rt, rp, erp, murp, err, murr, X_til_pt);


[S33, S34, S43, S44] = GSM(N, N, F(k), r(2), r(1), erp, murp, err, murr, X_til_rp);
[S11, S12, S21, S22] = GSM(N, N, F(k), r(3), r(2), ert, murt, erp, murp, X_til_pt);


Sl = SL(r(2), F(k), N, L);


I = eye(length(N), length(N));

U1 = inv(I - S22 * Sl * S33 * Sl);
U2 = inv(I - S33 * Sl * S22 * Sl);

STT(k, :, :) = (S11 + S12 * Sl * U2 * S33 * Sl * S21);
STR(k, :, :) = (S12 * Sl * U2 * S34);
SRT(k, :, :) = (S43 * Sl * U1 * S21);
SRR(k, :, :) = (S44 + S43 * Sl * U2 * S22 * Sl * S34);


end
%% Plots


% save('Stt4_ratio_1_modes_10_1mm', 'STT');
% save('Str4_ratio_1_modes_10_1mm', 'STR');
% save('Srt4_ratio_1_modes_10_1mm', 'SRT');
% save('Srr4_ratio_1_modes_10_1mm', 'SRR');

end