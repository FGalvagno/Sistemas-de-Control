% ITEM 1------------------------------------------------------------------------
%Generamos señal de entrada:
%Onda cuadrada alterna de 2ms de periodo y 12V de amplitud
pkg load control
[u,t] = gensig("square",2e-03,4e-03);
u=24*u-12;
R=47; L=1e-06; C=100e-09; %Parámetros del circuito RLC
A1 = [-R/L -1/L; 1/C 0];
B1 = [1/L; 0];
C1 = [1 0];
C2 = [0 1] 
D1 = [0];


% Creamos modelo de espacio de estados
sys = ss(A1,B1,C1,D1);
[i,t]=lsim(sys,u,t);

sys = ss(A1,B1,C2,D1);
[v_c,t] = lsim(sys,u,t);

figure('Name','1.1')
 subplot(3,1,1);plot(t,i);grid on; title('Corriente, I'); 
 subplot(3,1,2);plot(t,v_c);grid on; title('Tensión del capacitor, V_c');
 subplot(3,1,3);plot(t,u);grid on; title('Tensión de entrada, V_e');
 