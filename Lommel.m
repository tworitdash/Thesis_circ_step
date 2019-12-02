function [X] = Lommel(a, b, beta_rhoNu, beta_rhoMu, Nu, Mu)

    Xa = Zeta(a, beta_rhoNu, beta_rhoMu, Nu, Mu);
    Xb = Zeta(b, beta_rhoNu, beta_rhoMu, Nu, Mu);

    X = Xb - Xa;

end