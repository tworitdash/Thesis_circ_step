clear;

F = 4e9:0.01e9:21e9;

rp = 0.0405319403216/2;   % radius of the bigger waveguide
rr = 0.0405319403216/2.1; % radius of the smaller waveguide

Nr = 1:1:20;
Np = 1:1:20;

erp = 1;
murp = 1;
err = 1;
murr = 1;


[X_til_pr] = Inner_p(Nr, Np, rp, rr, erp, murp, err, murr);

parfor k = 1:length(F)
    disp('Freq iteration');
    disp(k);

[Spp_, Spr_, Srp_, Srr_] = GSM(Nr, Np, F(k), rp, rr, erp, murp, err, murr, X_til_pr);

    Spp(k, :, :) =  Spp_;
    Spr(k, :, :) =  Spr_;
    Srp(k, :, :) =  Srp_;
    Srr(k, :, :) =  Srr_;

end

data5 = read(rfdata.data,'S_Feko_5modes_each.s10p');
s_params_5 = extract(data5,'S_PARAMETERS');

figure;

plot(F * 1e-9, db(abs(squeeze(s_params_5(6, 6, :))))/2, 'LineWidth', 2); grid on;
hold on;
plot(F * 1e-9, db(abs(squeeze(Spp(:, 1, 1))))/2, 'LineWidth', 2); grid on;
