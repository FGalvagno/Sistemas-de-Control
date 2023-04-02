clc;clear all;close all;
% Generamos señal de entrada:
% Onda cuadrada alterna de 2ms de periodo y 12V de amplitud
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

%Cargamos los datos medidos
data = readtable('Curvas_Medidas_RLC.xls')
data = renamevars(data, ["Var1","Var2","Var3"], ["t","i","v_c"])

figure('Name', '1.2')
 subplot(3,1,1);plot(data.t,data.i);grid on; title('Corriente1');
 subplot(3,1,2);plot(data.t,data.v_c);grid on; title('Tensión en el capacitor, V_c1');

%Determinacion del modelo a partir de la ecuacion de 2do orden


