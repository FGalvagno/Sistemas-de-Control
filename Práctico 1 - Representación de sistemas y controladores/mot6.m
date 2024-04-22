% ITEM 6-------------------------------------------------------------------
clear all; close all; clc

Laa=5.28e-04;
J=2.61e-09;
Ra=20;
B=0;
Ki=0.0102;
Km=0.0605;

A = [-Ra/Laa -Km/Laa 0; Ki/J -B/J 0; 0 1 0];
B = [1/Laa 0; 0 -1/J; 0 0];
C = [1 0 0; 0 1 0; 0 0 1];
D = [0 0; 0 0; 0 0]
U = [12; 0]

X = [0; 0; 0];
Xp = [0; 0; 0;]
ii = 0;
t_etapa = 1e-7;
tF = 5e-3;
ref = 1;

%Constantes del PID
%Kp=0.1;Ki=0.01;Kd=5;color_='r';
Kp=300;Ki=100;Kd=0.01;color_='b';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;

N = tF/t_etapa
omega=zeros(1,N);
theta=zeros(1,N);
ia=zeros(1,N);
e=zeros(1,N);
acc=zeros(1,N);
ij = 0;

for ii= 1:N-1
    ij=ij+t_etapa;
    if (ij>=2.5e-3)
        U(2)=1e-3;
    end
    k = ii+2;
    Xp=A*X+B*U;
    X=X+Xp*t_etapa;
    Y=C*X+D*U;
    e(k) = ref-Y(3);
    U(1) = U(1)+A1*e(k)+B1*e(k-1)+C1*e(k-2);
    
    acc(ii+1)=U(1);
    ia(ii+1)=Y(1);
    omega(ii+1)=Y(2);
    theta(ii+1)=Y(3);
    
end

t=t_etapa:t_etapa:tF;
subplot(2,1,1);hold on;
plot(t,theta,color_);title('Salida: Posición, \theta_t');
subplot(2,1,2);hold on;
plot(t,omega,color_);title('Salida: Velocidad, \omega_t');
xlabel('Tiempo [Seg.]');


