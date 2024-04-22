clear all;close all;

X = -[0; 0; 0];
ii = 0;
t_etapa = 1e-7;
w_r = 1;
tF = 1e-3;
wRef = 1;
%Constantes del PID
Kp=.500;Ki=0.001;Kd=0.0001;color_='r';
% Kp=1;Ki=0;Kd=0.0001;color_='k';
% Kp=10;Ki=0;Kd=0;color_='b';
Ts=t_etapa;
A1=((2*Kp*Ts)+(Ki*(Ts^2))+(2*Kd))/(2*Ts);
B1=(-2*Kp*Ts+Ki*(Ts^2)-4*Kd)/(2*Ts);
C1=Kd/Ts;
e=zeros(tF/t_etapa,1);u=0;
for t=0:t_etapa:tF
ii=ii+1;k=ii+2;
X=modmotor(t_etapa, X, u);
e(k)=wRef-X(1); %ERROR
u=u+A1*e(k)+B1*e(k-1)+C1*e(k-2); %PID
x1(ii)=X(1);%Omega
x2(ii)=X(2);%wp
acc(ii)=u;
end
t=0:t_etapa:tF;
subplot(2,1,1);hold on;
plot(t,x1,color_);title('Salida y, \omega_t');
subplot(2,1,2);hold on;
plot(t,acc,color_);title('Entrada u_t, v_a');
xlabel('Tiempo [Seg.]');
% % Para verificar
% Laa=366e-6;
% J=5e-9;
% Ra=55.6;
% B=0;
% Ki=6.49e-3;
% Km=6.53e-3;
% num=[Ki]
% den=[Laa*J Ra*J+Laa*B Ra*B+Ki*Km ]; %wpp*Laa*J+wp*(Ra*J+Laa*B)+w*(Ra*B+Ki*Km)=Vq*Ki
% sys=tf(num,den)
% step(sys)