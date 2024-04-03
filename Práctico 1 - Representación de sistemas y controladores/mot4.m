% ITEM 4------------------------------------------------------------------------
%Asignamos los valores correspondientes al motor
clear all;
close all;
Laa=366e-6;
J=5e-9;
Ra=55.6;
B=0;
Ki=6.49e-3;
Km=6.53e-3;
delta=10e-07;
ts=5;

%Generamos la función de entrada de tensión Va y TL para el motor.

Va=12;
t=0:delta:(ts-delta);
T_L = t>=0;
for i=1:ts/delta-delta
 T_L(i)=0;
end

%Realizamos la integración por euler
ia = zeros(1,ts/delta);
wr = zeros(1,ts/delta);
o = zeros(1,ts/delta);
ia(1)=0; wr(1)=0; o(1)=0; %Cond. iniciales
Va=12;

for i=2:(ts/delta-1)
 ia(i)=ia(i-1)+delta*(-Ra*ia(i-1)/Laa-Km*wr(i-1)/Laa+Va/Laa);
 wr(i)=wr(i-1)+delta*(Ki/J*ia(i-1)-B/J*wr(i-1)-T_L(i-1)/J);
 o(i) = o(i-1) + delta*wr(i-1);
end

figure('Name','Torque vs Corriente')
 subplot(3,1,1);plot(t,wr);grid on; title('Velocidad angular, W_r'); 
 subplot(3,1,2);plot(t,ia);grid on; title('Corriente de armadura, I_a');

 %El torque máximo viene dado por la corriente máxima del motor (en el arranque)
 ia_max = max(ia)
 wr_max = max(wr)
 T_L_max = Ki*ia_max-B*wr_max
 
 %Se comprueba que son los valores máximos al aplicar un torque equivalente en
 %el modelo y observar que la velocidad angular se reduzca a cero.
 
 T_L = t>=0;
for i=1:ts/delta/2
 T_L(i)=0;
end
T_L = T_L*0.0014;
ia(1)=0; wr(1)=0; o(1)=0; %Cond. iniciales
Va=12;
for i=2:(ts/delta-1)
 ia(i)=ia(i-1)+delta*(-Ra*ia(i-1)/Laa-Km*wr(i-1)/Laa+Va/Laa);
 wr(i)=wr(i-1)+delta*(Ki/J*ia(i-1)-B/J*wr(i-1)-T_L(i-1)/J);
 o(i) = o(i-1) + delta*wr(i-1);
end

figure('Name', 'W_r al aplicar T_L máximo')
plot(t, wr)
 
 
