%% Calculo de parametros 2do orden
% Especificaciones de diseño
S = 15;
ts = 4;
e = 0;
Tm=0.19

psita = -log(S/100)/sqrt(pi^2+log(S/100)^2);
wo = 4/(psita*ts);
wd = wo*sqrt(1-psita^2);
td = 2*pi/wd;
m = round(td/Tm);

r = exp(-psita*wo*Tm)
Omega = wd*Tm

[x y] = pol2cart(Omega,r);
z1 = x+j*y;
z2 = x-j*y;

G=zpk(-10,[-3 0],10)
Gd=c2d(G,Tm,'zoh')

sisotool(Gd)

C
F=feedback(C*Gd,1) % sistema de lazo cerrado
pole(F)
zero(F)
pzmap(F)
step(F) % respuesta al escalon

%parametros PID
syms kp kd ki k c b 
eq1 = kp + kd + ki == k
eq2 = (kp+2*kd)/k == b
eq3 = kd/k == c
S1=solve(eq1,eq2,eq3,kp,ki,kd);

c_o = zero(C)
p_o = pole(C)
k_o = C.K

b = c_o(1)+c_o(2)
c = c_o(1)*c_o(2) 
k = C.K

kp = b*k - 2*c*k
kd = c*k
%ki = k*(c - b + 1)
ki = 0
