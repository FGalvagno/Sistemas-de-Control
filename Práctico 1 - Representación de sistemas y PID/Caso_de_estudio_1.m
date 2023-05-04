%1------------------------------------------------------------------------
%Generamos señal de entrada:
%Onda cuadrada alterna de 2ms de periodo y 12V de amplitud
[u,t] = gensig("square",2e-03,4e-03);
u=24*u-12;
R=4.7e03; L=10e-06; C=100e-09; %Parámetros del circuito RLC
A1 = [-R/L -1/L; 1/C 0];
B1 = [1/L; 0];
C1 = [R 0];
D1 = [0];


% Creamos modelo de espacio de estados
sys = ss(A1,B1,C1,D1);

[i,t]=lsim(sys,u,t);

figure('Name','1.1')
 subplot(2,1,1);plot(t,i);grid on; title('Corriente, I'); 
 subplot(2,1,2);plot(t,u);grid on; title('Tensión de entrada, V_e');

%2------------------------------------------------------------------------
% Modelamos de nuevo el mismo sistema, pero con valores de R diferentes
R = 5.6e03;

A2 = [-R/L -1/L; 1/C 0];
B2 = [1/L; 0];
C2 = [R 0];
D2 = [0];


% Creamos modelo de espacio de estados
sys = ss(A2,B2,C2,D2);

[i,t]=lsim(sys,u,t);

figure('Name', '1.2')
 subplot(2,1,1);plot(t,i);grid on; title('Corriente, I');
 subplot(2,1,2);plot(t,u);grid on; title('Tensión de entrada, V_e');
 
%3------------------------------------------------------------------------
clear all;
%Cargamos los datos medidos
datos = xlsread('Curvas_Medidas_RLC.xls');
t = datos(:,1);
i = datos(:,2);
v_c = datos(:,3);

figure('Name', '1.2')
 subplot(2,1,1);plot(t,i);grid on; title('Corriente1');
 subplot(2,1,2);plot(t,v_c);grid on; title('Tensión en el capacitor, V_c1');

% Determinamos el modelo del circuito RLC a partir del método desarollado
% por Chen

td = 0.01;
t1 = 0.002;
k=max(v_c);
[val lugar] = min(abs((t1+td)-t))

k1 = v_c(lugar)/k-1
[val lugar] = min(abs((2*t1+td)-t))
k2 = v_c(lugar)/k-1
[val lugar] = min(abs((3*t1+td)-t))
k3 = v_c(lugar)/k-1

b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1 = (k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2 = (k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta= (2*k1^3+3*k1*k2+k3-sqrt(b))/sqrt(b);

T1 =-t1/log(alfa1);
T2 =-t1/log(alfa2);
%T3 = beta*(T1-T2)+T1

s = tf('s');

%sys_id = k*(T3*s+1)/(T1*s+1)/(T2*s+1)
sys_id = k/(T1*s+1)/(T2*s+1)

figure
    step(sys_id)
    
[N D] = tfdata(sys_id, "v")
L = 0.1; %Definimos un valor para L
C = D(1)/L;
R = D(2)/C;

%4------------------------------------------------------------------------
[u,t] = gensig("square",0.08,0.16);
u = -(24*u-12);
A4 = [-R/L -1/L; 1/C 0];
B4 = [1/L; 0];
C4 = [1, 0];
D4 = [0];

% Creamos modelo de espacio de estados
sys_est = ss(A4,B4,C4,D4);

[i_est,t_sim]=lsim(sys_est,u,t);
 
figure('Name','1.4')
hold on
plot(t, i_est)
plot(datos(:,1)-0.01,datos(:,2));grid on; title('Corriente');
legend('RLC estimado','RLC medido')
xlim([0 0.08])