function G = segmap(s,S)

% s = Gauss Quadrature nodes
% S = triangle nodes
psi1 = 0.5*(1 - s);
psi2 = 0.5*(1 + s);

w11 = S(1,1);
w12 = S(1,2);

w21 = S(2,1);
w22 = S(2,2);

G1 = w11*psi1 + w21*psi2;
G2 = w12*psi1 + w22*psi2;

G = [G1 G2];
