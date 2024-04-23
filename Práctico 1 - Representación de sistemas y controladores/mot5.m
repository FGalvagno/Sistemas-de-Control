% ITEM 5------------------------------------------------------------------------
clear all; close all;
% Cargamos los datos medidos
datos = xlsread('Curvas_Medidas_Motor_2024.xlsx');
t = datos(:,1);
w_r= datos(:,2);
i_a = datos(:,3);
u = datos(:,4);
t_l = datos(:,5);

figure('Name', 'Ítem 5')
subplot(4,1,1);plot(t,i_a);grid on; title('Corriente de Armadura, i_a');
subplot(4,1,2);plot(t,w_r);grid on; title('Velocidad Angular, W_r');
subplot(4,1,3);plot(t,u);grid on; title('Tensión de entrada, V_i');
subplot(4,1,4);plot(t,t_l);grid on; title('Torque de carga, T_L');

%% Obtención de la funcion de transferencia del motor
syms V Ra La I Ki w Jm Km Bm Tl s real
eq1=s*I==-Ra/La*I-Km/La*w+1/La*V;
eq2=s*w==Ki/Jm*I-Bm/Jm*w;
S1=solve(eq1,eq2,w,V);
wr_va=collect(S1.w/S1.V,s);


syms V Ra La I Ki w Jm Km Bm Tl s real
eq1=s*I==-Ra/La*I-Km/La*w;
eq2=s*w==Ki/Jm*I-Bm/Jm*w-Tl/Jm;
S1=solve(eq1,eq2,w,Tl);
wr_tl=collect(S1.w/S1.Tl,s);


% tomamos Bm = 0
wr_va = subs(wr_va, Bm, 0)
wr_tl = subs(wr_tl, Bm, 0)
pretty(wr_va);
pretty(wr_tl);
% Según TVF tenemos
limit(wr_va, s, 0)
limit(wr_tl, s, 0)
%% Determinamos las funciones de transferencia a partir del método desarollado
% por Chen

% td: tiempo de delay de entrada
% t1: tiempo inicial para el algoritmo

td = 0.0351;
t1 = 0.0001;
%k=max(w_r);
k = 198
[val lugar] = min(abs((t1+td)-t))

 k1 = w_r(lugar)/k-1;
 [val lugar] = min(abs((2*t1+td)-t));
 k2 = w_r(lugar)/k-1;
 [val lugar] = min(abs((3*t1+td)-t));
 k3 = w_r(lugar)/k-1;
%k1 = 135.58/k-1
%k2 = 191.60/k-1
%k3 = 198.30/k-1

b=4*k1^3*k3-3*k1^2*k2^2-4*k2^3+k3^2+6*k1*k2*k3;
alfa1 = (k1*k2+k3-sqrt(b))/(2*(k1^2+k2));
alfa2 = (k1*k2+k3+sqrt(b))/(2*(k1^2+k2));
beta= (2*k1^3+3*k1*k2+k3-sqrt(b))/sqrt(b);

T1 =-t1/log(alfa1);
T2 =-t1/log(alfa2);
%T3 = beta*(T1-T2)+T1
s = tf('s');
G1=(k/12)/(T1*s +1)/(T2*s +1); %G1 normalizada.


figure
hold on
step(G1*12,5e-4, 'r') %Simulo aplicando 12V a la entrada.
plot(datos(:,1)-0.0351,datos(:,2));grid on; title('Velocidad Angular W_r');
legend('FdT estimado','FdT medido')
hold off


%% Obtención de parámetros del motor
%Aproximación de la resistencia de armadura 28.13 en una primera aprox.
%Bajamos la resistencia en otra aproximación a 20 Ohm según el primer pico
%de la corriente
%Ra = max(u)/max(i_a)
Ra = 20
Km = max(u)/max(w_r)
%A partir de la ganancia en estado estable de Wr/Tl despejamos Ki 
K_ss = (198-164)/1.05e-3
Ki = Ra/Km/K_ss
[N D] = tfdata(G1, 'v')

%Normalizamos el numerador y el denominador
fac = Ki/N(3)
N = N*fac
D = D*fac

B = 0; 
J = D(2)/Ra
Laa = D(1)/J


%%
%Verificamos el modelo aproximado vs las mediciones
delta=1e-6;
ts=0.6;

%Generamos la función de entrada de tensión Va y TL para el motor.
t_s =0:delta:(ts-delta);
ia = zeros(1,ts/delta);
wr = zeros(1,ts/delta);
o = zeros(1,ts/delta);
va = zeros(1,ts/delta);
T_L = zeros(1,ts/delta);
ia(1)=0; wr(1)=0; o(1)=0; %Cond. iniciales

for i=1:(ts/delta)
    if i*delta < 0.0351
        va(i) = 0;
    else
        va(i) = 12;
    end
end


 for i=1:(ts/delta)
     if i*delta < 0.1863
         T_L(i) = 0;
     elseif i*delta < 0.3372
         T_L(i) = 1.03e-3;
     elseif i*delta < 0.4866
         T_L(i) = 0;
     else
         T_L(i) = 1.03e-3;
     end
 end
        
for i=2:(ts/delta-1)
 ia(i)=ia(i-1)+delta*(-Ra*ia(i-1)/Laa-Km*wr(i-1)/Laa+va(i-1)/Laa);
 wr(i)=wr(i-1)+delta*(Ki/J*ia(i-1)-B/J*wr(i-1)-T_L(i-1)/J);
 o(i) = o(i-1) + delta*wr(i-1);
end

figure('Name','Torque vs Corriente')
 subplot(3,1,1);
 plot(t_s,wr);grid on; title('Velocidad angular, W_r'); 
 hold on
 plot(t,w_r);
 legend('Aprox', 'Medido')
 hold off
 
 subplot(3,1,2);plot(t_s,ia);grid on; title('Corriente de armadura, I_a');
 hold on
 plot(t,i_a);
 hold off
 
 
 subplot(3,1,3);plot(t_s,T_L);grid on; title('Torque de carga, T_L');
 hold on
 plot(t,t_l);
 hold off
 
 