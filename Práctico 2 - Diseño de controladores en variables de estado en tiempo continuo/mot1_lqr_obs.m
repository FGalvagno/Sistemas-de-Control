clear all; close all; clc
Laa=5.28e-04;
J=2.61e-09;
Ra=20;
Bm=0;
Ki=0.0102;
Km=0.0605;
%% Modelo del motor en espacio de estados

A = [-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B = [1/Laa; 0; 0];
C = [0 1 0; 0 0 1];
D = [0 0; 0 0; 0 0];

Q = diag([62500 2.5e-5  1  10])
R = 1;

% Matriz ampliada 
Aa = [A zeros(3,1); -C(2,:) 0]
Ba = [B(:,1); 0]

%Q=diag([100 50 10 ]);R=1e-1;
%K = lqr(Aa,Ba,Q,R)

 %Construcción del Hamiltoniano para el cálculo del controlador
Ha=[Aa -Ba*inv(R)*Ba'; -Q -Aa'];
[n,va]=size(Ha);
[V,D]=eig(Ha);MX1X2=[];
for ii=1:n
 if real(D(ii,ii))<0
 MX1X2=[MX1X2 V(:,ii)];
end
end
MX1=MX1X2(1:n/2,:); MX2=MX1X2(n/2+1:end,:);
 % P=abs(MX2*inv(MX1));
P=real(MX2*inv(MX1));
K=inv(R)*Ba'*P;
disp('Controlador ampliado en ')
eig(Aa-Ba*K)
%% Calculo del observador
Ao = A';
Bo = C';
Co = B';
Qo = diag([1e1 10 .09]);
Ro = diag([1e1 1e1]);
Ko = lqr(Ao,Bo,Qo,Ro)';
%poles_o=eig(A-Ko*C(2:3,:))
disp('Observador en ')
eig(A-Ko*C)
%% Calculo del tiempo de simulacion
poles = eig(Aa-Ba(:,1)*K)
lambda = max(poles);
tr = abs(log(0.95)/lambda);
delta = 2.5e-7;
ts = 0.5;

%% Señales de entrada
t_s =0:delta:(ts-delta);
T_L = zeros(1,ts/delta);
ia = zeros(1,ts/delta);
u = zeros(1,ts/delta);
wr = zeros(1,ts/delta);
o_ref = zeros(1,ts/delta);
o = zeros(1,ts/delta);

ia_o = zeros(1,ts/delta);
wr_o = zeros(1,ts/delta);
o_o = zeros(1,ts/delta);
e_o = zeros(2,ts/delta);
va = zeros(1,ts/delta);
ia(1)=0; wr(1)=0; o(1)=0; %Cond. iniciales
X_hat = [ia(1) wr(1) o(1)]';
x_hat = [0 0 0]'; 

 for i=1:(ts/delta)
     if i*delta < 5
         o_ref(i) = -pi/2;
     elseif i*delta < 10
         o_ref(i) = pi/2;
     elseif i*delta < 15
         o_ref(i) =-pi/2;
     else
         o_ref(i) = pi/2;
     end
 end
 

T_L = 1e-3;
%% Simulación

ia(1) = 0;
o(1) = 0;
wr(1) = 0;
%T_L =  1.03e-3;
%T_L =  0;
X = [ia(1); wr(1); o(1)];
u(1) = 0;
index = 0;
int = 0;
zeta = 0;
for i = 1:(ts/delta-1)
      zeta_p = o_ref(i)-C(2,:)*X;
      zeta = zeta + zeta_p*delta;
      u(i) = -K*[X;zeta];
      x1_p = -Ra*X(1)/Laa-Km*X(2)/Laa+u(i)/Laa;
      x2_p = Ki*X(1)/J-Bm*X(2)/J-T_L/J;
      x3_p = X(2);
      xP = [x1_p x2_p x3_p]';
      X = X+delta*xP;
      
      ia(i) = X(1);
      wr(i) = X(2);
      o(i) = X(3);
      X = [ia(i) wr(i) o(i)]';

end
x_hat=[0;0;0]
% figure
% hold on
% plot(t_s, o);
% plot(t_s, o_ref, 'Color', 'g')
% title('Salida del sistema \theta_r(t)')
% xlabel('Tiempo [s]')
% ylabel('Ángulo [rad]')
% grid
% hold off
X = [0;0;0]

zeta = 0;
y = 0;
y_o = 0;
for i = 1:(ts/delta)    
      y = C*X;
      zeta_p = o_ref(i)-C(2,:)*X;
      zeta = zeta + zeta_p*delta;
      u(i) = -K*[x_hat;zeta];
      
      x1_p = -Ra*X(1)/Laa-Km*X(2)/Laa+u(i)/Laa;
      x2_p = Ki*X(1)/J-Bm*X(2)/J-T_L/J;
      x3_p = X(2);
      xP = [x1_p x2_p x3_p]';
      X = X+delta*xP;
      
      y_o = C*x_hat;
      x_hatp = A*x_hat+B*u(i)+Ko*(y-y_o);
      x_hat = x_hat + delta*x_hatp;
      
      ia_o(i) = x_hat(1);
      wr_o(i) = x_hat(2);
      o_o(i) = x_hat(3);
      e_o(:,i) = y-y_o;

end

figure
subplot(3,1,1)
hold on
plot(t_s, o);
plot(t_s, o_ref, 'Color', 'g')
plot(t_s, o_o, 'Color', 'y')
title('Salida del sistema \theta_r(t)')
xlabel('Tiempo [s]')
ylabel('Ángulo [rad]')
grid
hold off


subplot(3,1,2)
plot(t_s,u,'LineWidth',1.5);
title('Acción de control')
xlabel('Tiempo [s]')
grid

subplot(3,1,3)
hold on
plot(t_s,ia,'LineWidth',1.5);
title('Salida del sistema i_a(t)')
xlabel('Tiempo [s]')
ylabel('Corriente [A]')
grid
hold off

figure
hold on
plot(t_s, e_o(1,:))
plot(t_s, e_o(2,:), 'Color', 'g')
title('Error en la observacion')