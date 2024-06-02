%% Tarea 2
% *De la especificaciones de diseño dadas archivo adjunto (sobrepaso, tiempo 
% de establecimiento 2 % y error de régimen ante escalón=0)*

%% Calculo de parametros 2do orden
% Especificaciones de diseño
S = 15;
ts = 4;
e = 0;
Tm=0.19;

psita = -log(S/100)/sqrt(pi^2+log(S/100)^2)
wo = 4/(psita*ts)
wd = wo*sqrt(1-psita^2)
td = 2*pi/wd
m = round(td/Tm);

r = exp(-psita*wo*Tm)
Omega = wd*Tm

[x y] = pol2cart(Omega,r);
z1 = x+j*y
z2 = x-j*y

c = -10;
p1 = -3;
p2 = 0;
K = 10;

G=zpk(-10,[-3 0],10)
Gd=c2d(G,Tm,'zoh')

%sisotool(Gd);
%% Diseño de un controlador PI
% Se diseña con _sisotool_ un controlador PI para los parámetros de diseño. 
%%
load("pi_out.mat");
C
%% 
% Simulando el controlador obtenido

F=feedback(C*Gd,1) % sistema de lazo cerrado
pzmap(C*Gd)
pole(F)
zero(F)
pzmap(F)
step(F) % respuesta al escalon
stepinfo(F)
%% 
% Se puede observar que a partir del diseño realizado, no se cumple el sobrepaso 
% requerido. Si cancelamos el polo de la planta con el cero del controlador PI, 
% no se pueden cumplir las especificaciones de diseño a partir de un cambio de 
% ganancia del controlador.
%% Diseño de un controlador PD
% Puesto que el sistema posee un integrador puro, su error en estado estable 
% para una entrada escalón será nulo (sistema tipo 1). Se decide diseñar un control 
% PD para mejorar su respuesta transitoria.
%%
%sisotool(Gd)
load('pd_out.mat')
F=feedback(C*Gd,1) % sistema de lazo cerrado
pole(F)
zero(F)
pzmap(F)
step(F) % respuesta al escalon
stepinfo(F)
%% 
% Se puede observar que a partir del diseño del controlador PD se cumplen 
% las especificaciones de diseño y además se obtiene sobrepaso nulo.
%% Simulacion del PD mediante SIMULINK
% Para simular el controlador PD, se deben de calcular los parámetros del controlador 
% KD y KP.
%%
%parametros PID
syms kp kd ki k c b 
eq1 = kp + kd == k
eq2 = kd/(kp+kd) == b
S1=solve(eq1,eq2,kp,kd);

c_o = zero(C);
b = c_o(1);
c = 0;
k = C.K;

kp = k - b*k
kd = b*k
ki= 0;
T = Tm;
sim('PID_digital_tarea.slx')


subplot(5,1,1);
plot(tout,yout(:,1)) 
grid on
title('Salida')
subplot(5,1,2);
plot(tout,yout(:,2)) 
grid on
title('Señal de control')
subplot(5,1,3);
plot(tout,yout(:,3)) 
grid on
title('Error')
subplot(5,1,4);
plot(tout,yout(:,4)) % salida del sistema
grid on
title('Acción derivativa')
subplot(5,1,5);
plot(tout,yout(:,6)) % salida del sistema
grid on
title('Acción proporcional')