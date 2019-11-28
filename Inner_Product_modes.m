
%% Zeros of Bessel's function (For TM) and zeros of the derivative of the Bessel's function (For TE) Calculation
clear;
% [m, n, xm] = meshgrid(1:1:50, 1:1:50, linspace(0.1, 10000, 100000));
m = 1:1:50;
mode = "TE";

for j = 1:length(m)

xm = linspace(0.1, 10000, 10000);

Jm = @(z) besselj(m(j), z);                                         % Bessel's function

dz = 1e-5;

Jm_der = @(z) (besselj(m(j), z + dz) - besselj(m(j), z - dz))./(2 * dz); % Derivative of the Bessel's function

if mode == "TE"
    fun = Jm_der;
else
    fun = Jm;
end

ym = fun(xm);                   % Values of the Bessel's function at the positions defined by xm

chsign = find(diff(sign(ym)));  % Detection of the sign changes 

% xmn_all = fzero(fun, xm(chsign));

% xmn_all = zeros(size(chsign));

for i = 1:size(chsign, 2)
    xmn_all(j, i) = fzero(fun, xm(chsign(i)));  % finding the roots near the points where Jm changes sign
end


end
%xmn = xmn_all(n); 


% beta_rho = xmn/r;  % wave number along the rho direction (\beta_{\rho})

%% Inner Product Calculation

X_til = zeros(Np(end), Nr(end));

for p = 1:length(Np)
    for r = 1:length(Nr)
        disp(p);
        disp(r);
        grad_Phi_rhop = Nup(p) .* cos(mp .* phir_) .* besselj_der(mp, beta_rhop(p) .* rhor_) .* beta_rhop(p);
        grad_Phi_phip = (-1./rhor_) .* Nup(p) .* mp .* sin(mp .* phir_) .* besselj(mp,  beta_rhop(p) .* rhor_);
        
        grad_Phi_rhor = Nur(r) .* cos(mr .* phir_) .* besselj_der(mr, beta_rhor(r).* rhor_) .* beta_rhor(r);
        grad_Phi_phir = (-1./rhor_) .* Nur(r) .* mr .* sin(mr .* rhor_) .* besselj(mp,  beta_rhor(r) .* rhor_);
        
        if (modep == "TE" && moder == "TE") || (modep == "TM" && moder == "TM")
            X_til_pr = (grad_Phi_rhop .* grad_Phi_rhor +  grad_Phi_phip .* grad_Phi_phir)...
                .* rhor_ .* drho .* dphi;
            X_til(p, r) = sum(sum(X_til_pr));
        elseif (modep == "TE" && moder == "TM")
            X_til(p, r) = 0;
        elseif (modep == "TM" && moder == "TE")
            X_til_pr = (grad_Phi_rhop .* grad_Phi_phir - grad_Phi_rhor .* grad_Phi_rhop)...
                .* rhor_ .* drho .* dphi;
            X_til(p, r) = sum(sum(X_til_pr));
        else
        end
    end
end
