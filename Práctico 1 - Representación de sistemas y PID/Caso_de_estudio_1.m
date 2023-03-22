clc;clear all;
% Generamos señal de entrada:
% Onda cuadrada alterna de 2ms de periodo y 12V de amplitud
[u,t] = gensig("square",2e-03,4e-03);
u=24*u-12;
R=4.7e03; L=10e-06; C=100e-09; %Parámetros del circuito RLC
A = [-R/L -1/L; 1/C 0]
B = [1/L; 0]
C= [R 0]
D= [0]


% Creamos modelo de espacio de estados
sys = ss(A,B,C,D)


 [i,vc]=lsim(A,B,C,D,u,t)
 
 subplot(3,2,1);plot(t,i);grid on; title('Corriente'); hold on
 subplot(3,2,2);plot(t,vc);grid on; title('Tensión en el capacitor, V_c'); hold on
 subplot(3,2,3);plot(t,u);grid on; title('Tensión de entrada, V_e'); hold on
 subplot(3,1,1);plot(t,i);grid on; title('Corriente'); hold on
 subplot(3,2,2);plot(t,vc);grid on; title('Tensión en el capacitor, V_c'); hold on
 subplot(3,2,3);plot(t,u);grid on; title('Tensión de entrada, V_e'); hold on