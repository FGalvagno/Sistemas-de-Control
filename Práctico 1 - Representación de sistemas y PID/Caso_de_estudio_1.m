clc;clear all;close all;
% Generamos señal de entrada:
% Onda cuadrada alterna de 2ms de periodo y 12V de amplitud
%{
[u,t] = gensig("square",2e-03,4e-03);
u=24*u-12;
R=4.7e03; L=10e-06; C=100e-09; %Parámetros del circuito RLC
A1 = [-R/L -1/L; 1/C 0];
B1 = [1/L; 0];
C1 = [R 0];
D1 = [0];


% Creamos modelo de espacio de estados
sys = ss(A1,B1,C1,D1);

[i,vc]=lsim(A1,B1,C1,D1,u,t);

figure('Name','1.1')
 subplot(3,1,1);plot(t,i);grid on; title('Corriente');
 subplot(3,1,2);plot(t,vc);grid on; title('Tensión en el capacitor, V_c'); 
 subplot(3,1,3);plot(t,u);grid on; title('Tensión de entrada, V_e');

% Modelamos de nuevo el mismo sistema, pero con valores de R diferentes
R = 5.6e03;

A2 = [-R/L -1/L; 1/C 0];
B2 = [1/L; 0];
C2 = [R 0];
D2 = [0];


% Creamos modelo de espacio de estados
sys = ss(A2,B2,C2,D2);

[i,vc]=lsim(A2,B2,C2,D2,u,t);

figure('Name', '1.2')
 subplot(3,1,1);plot(t,i);grid on; title('Corriente1');
 subplot(3,1,2);plot(t,vc);grid on; title('Tensión en el capacitor, V_c1');
 subplot(3,1,3);plot(t,u);grid on; title('Tensión de entrada, V_e1');
%}
%Cargamos los datos medidos


datos = table2array(readtable('Curvas_Medidas_RLC.xls'));
t = datos(:,1);
v_c = datos(:,3);



%figure('Name', '1.2')
 subplot(3,1,1);plot(datos(:,1),datos(:,2));grid on; title('Corriente1');
 subplot(3,1,2);plot(datos(:,1),datos(:,3));grid on; title('Tensión en el capacitor, V_c1');

% Determinamos el modelo del circuito RLC a partir del metodo desarollado
% por Chen

td = 0.01
t1 = 0.002

k=max(v_c)

[val lugar] = min(abs((t1+td)-t));
k1 = v_c(lugar)/k-1
[val lugar] = min(abs((2*t1+td)-t));
k2 = v_c(lugar)/k-1
[val lugar] = min(abs((3*t1+td)-t));
k3 = v_c(lugar)/k-1

b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1 = (k1*k2+k3-sqrt(b))/(2*(k1^2+k2))
alfa2 = (k1*k2+k3+sqrt(b))/(2*(k1^2+k2))
beta= (2*k1^3+3*k1*k2+k3-sqrt(b))/sqrt(b)

T1 =-t1/log(alfa1)
T2 =-t1/log(alfa2)
%T3 = beta*(T1-T2)+T1


s = tf('s')

%sys_id = k*(T3*s+1)/(T1*s+1)/(T2*s+1)
sys_id = k*(T1*s+1)/(T2*s+1)
%V_c/V = 1/(LCs^2+R1*C*s+1)

step(sys_id)

