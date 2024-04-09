% ITEM 2------------------------------------------------------------------------
clear all; close all;
%Cargamos los datos medidos
datos = xlsread('Curvas_Medidas_RLC_2024.xlsx');
t = datos(:,1);
i = datos(:,2);
v_c = datos(:,3);

figure('Name', '1.2')
 subplot(2,1,1);plot(t,i);grid on; title('Corriente1');
 subplot(2,1,2);plot(t,v_c);grid on; title('Tensión en el capacitor, V_c1');

%Del gráfico aproximamos R a partir de Vmax/Imax, a partir del primer pico
%Raprox = 335 ohm
%Luego de contrastar los resultados, ajustamos los valores de R para llegar
%a una aproximación mejor de R = 270 ohm
R = 270
 
 
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
    
[N D] = tfdata(sys_id, 'v')

C = D(2)/R
L = D(1)/C
