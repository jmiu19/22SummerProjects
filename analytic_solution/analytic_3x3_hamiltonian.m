syms E G C X R x;

M=[E+1i*G,	             C,        R;	         
        C,          E-1i*G,        R;    
        R,	             R,        X];
xI=[x, 0, 0;
    0, x, 0;
    0, 0, x];    

characEq = det(M-xI);
eigVals = eig(M);
%eigValsLaTeX = latex(eigVals)



<<<<<<< HEAD
% C = 8;
% detuning = 7.5;
% R = 5;
% G = 4.180513;
% G = 7;
% X = 1650;
% E = X + detuning;
=======
C = 8;
R = 10;
G = 5
E = 1650;
X = 1645;
>>>>>>> 501c2673fea82708740fd2746e76e0485560dec2


kappa = ((C*C)-(E*E)-(G*G)+(2*R*R)-(2*E*X))/(3) + ((2*E+X)^2)/9
sigma = C*R*R-E*R*R+((E*E*X+G*G*X-C*C*X)/2)-((2*E+X)*(-C*C+E*E+2*E*X+G*G-2*R*R))/(6)+((2*E+X)^3)/(27)
zeta = (sigma + sqrt(sigma*sigma - kappa^3))^(1/3)

mu1 = (kappa/zeta - zeta)
mu2 = (kappa/zeta + zeta)

val1 = (2*E+X)/3 + (kappa/zeta + zeta);
val2 = (2*E+X)/3 - (1/2)*(kappa/zeta + zeta) - (sqrt(3)/2)*(kappa/zeta - zeta)*1i;
val3 = (2*E+X)/3 - (1/2)*(kappa/zeta + zeta) + (sqrt(3)/2)*(kappa/zeta - zeta)*1i;


% check = [eval(eigVals(1))- val1, eval(eigVals(2))- val2, eval(eigVals(3))- val3];

A = C^2 - E^2 - G^2 + 2*R^2 - 2*E*X;
B = 2*E + X
D = C*R^2 - E*R^2 + ((E^2)*X + (G^2)*X - (C^2)*X)/2

diff = D^2 + (A*B*D)/3 + (A^2*B^2)/36 + (2*C*B^3)/27 - A^3/27 - (A^2*B^2)/81 -(2*A^2*B^2)/81





M2=[E+1i*G,	              C;	         
         C,          E-1i*G];
xI2=[x, 0;
     0, x];    

characEq2 = det(M2-xI2);
eigVals2 = eig(M2);
%eigVals2LaTeX = latex(eigVals2);