   c0 = 3e8;

    Str = load('Xmn.mat');
    Xmn = Str.Xmn;
    
    rr = 0.0405319403216/1.9;
    F = 1e9:0.5e9:21e9;
    
    beta_rho = (Xmn(1).xmn)./rr;
    
    beta = (2 * pi * F) ./ c0;
    
%     beta_z = -1j .* sqrt(-(beta.^2 - beta_rho.^2));
    beta_z = sqrt((beta.^2 - beta_rho.^2));
    
    L = 0.001;
    
    Phase = exp(-1j * beta_z * L);
    
    plot(F * 1e-9, angle(Phase)*180/pi, 'LineWidth', 2); grid on;
   