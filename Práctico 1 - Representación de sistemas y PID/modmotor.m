function [X]=modmotor(t_etapa, xant, accion)
% modmotor modela un motor de continua con Euler
% Entradas:
%   t_etapa: tiempo de integración
%   xant: X de la iteración anterior
%   accion: accion de control aplicada al motor
% Salidas:
%   X: un vector de dos elementos omega y wp
 Laa=366e-6; J=5e-9;Ra=55.6;B=0;Ki=6.49e-3;Km=6.53e-3; 
 Va=accion;
 h=1e-7;
 omega= xant(1);
 wp= xant(2);
for ii=1:t_etapa/h
  wpp =(-wp*(Ra*J+Laa*B)-omega*(Ra*B+Ki*Km)+Va*Ki)/(J*Laa);
  wp=wp+h*wpp;
  omega = omega + h*wp;
 end
X=[omega,wp];
end