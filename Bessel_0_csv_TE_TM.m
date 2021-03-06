%% Zeros of Bessel's function (For TM) and zeros of the derivative of the Bessel's function (For TE) Calculation
clear;
% [m, n, xm] = meshgrid(1:1:50, 1:1:50, linspace(0.1, 10000, 100000));
m = 1:1:1000;
mode = "TE";
%mode = "TM";

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
csvwrite('TE_Xmn', xmn_all); % use it when mode is TE
% csvwrite('TM_Xmn', xmn_all);   % use it when mode is TM


% beta_rho = xmn/r;  % wave number along the rho direction (\beta_{\rho})