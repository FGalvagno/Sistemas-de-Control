clear all; close all; clc
Laa=5.28e-04;
J=2.61e-09;
Ra=20;
Bm=0;
Ki=0.0102;
Km=0.0605;
%% Modelo del motor en espacio de estados

A = [-Ra/Laa -Km/Laa 0; Ki/J -Bm/J 0; 0 1 0];
B = [1/Laa 0; 0 -1/J; 0 0];
C = [1 0 0; 0 1 0; 0 0 1];
D = [0 0; 0 0; 0 0];

Q = diag([62500 2.5e-5  1  10])
R = 1;

% Matriz ampliada (referencia diferente de 0)
Aa = [A zeros(3,1); -C(3,:) 0]
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

%% Calculo del observador
Ao = A';
Bo = C(2:3,:)';
Co = B';
Qo = diag([1 10 100]);
Ro = diag([1 1]);
Ko = lqr(Ao,Bo,Qo,Ro);
Ko = Ko';
poles_o=eig(A-Ko*C(2:3,:))
%% Calculo del tiempo de simulacion
poles = eig(Aa-Ba(:,1)*K)
lambda = max(poles);
tr = abs(log(0.95)/lambda);
delta = 2.5e-7;
ts = 20;

%% Señales de entrada
t_s =0:delta:(ts-delta);
T_L = zeros(1,ts/delta);
ia = zeros(1,ts/delta);
u = zeros(1,ts/delta);
wr = zeros(1,ts/delta);
o_ref = zeros(1,ts/delta);
o = zeros(1,ts/delta);
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
x = [ia(1); wr(1); o(1)];
X = [ia(1);wr(1);o(1)];
u(2,:) = 0;
y(2,:) = 0;
Y_o(2,:) = 0;
index = 0;
int = 0;
for i = 1:(ts/delta-1)
      y(:,i) = C(2:3,:)*X;
      psi_p = o_ref(i)-C(3,:)*X;
      psi(i) = int + psi_p*delta;
      u(i) = -K(1:3)*X_hat-K(4)*psi(i);
      ia(i) = x(1);
      wr(i) = x(2);
      o(i) = x(3);
      x1_p = -Ra*x(1)/Laa-Km*x(2)/Laa+u(i)/Laa;
      x2_p = Ki*x(1)/J-Bm*x(2)/J-T_L/J;
      x3_p = x(2);
      xP = [x1_p x2_p x3_p]';
      x = x+delta*xP;
      Y_o(:,i) = C(2:3,:)*X_hat;
      
      x_hatp =A*x_hat+B*u(i)+Ko*(y(:,i)-Y_o(:,i));
      x_hat = x_hat + delta*x_hatp;
      
      int = psi(i);
      X = [ia(i) wr(i) o(i)]';
      
      ia_o(i) = x_hat(1);
      wr_o(i) = x_hat(2);
      o_o(i) = x_hat(3);
      X_hat = [ia_o(i) wr_o(i) o_o(i)]';
%     psi_p = o_ref(i)-C(3,:)*X;
%     psi(i) = int + psi_p*delta;
%     u(i) = -K(1:3)*X-K(4)*psi(i);
%     ia(i) = x(1);
%     wr(i) = x(2);
%     o(i) = x(3);
%     x = x+delta*xP;
%     
%     ia_hat(i)= xObs(1);
%     omega_hat(i)= xObs(2);
%     theta_hat(i)= xObs(3);
%     X = [ia(i) wr(i) o(i)]';
%     y(i) = 
%     y_hat(i) = C*x_hat
%     
%     X_hat = A*x_hat+B*u(i)+Ko*(y(:,i)-y_hat(:,i));
%     x_hat = x_hat + X_hat*h
%     int = psi(i)  ;
end

ia_o(i+1)=ia_o(i);
o_o(i+1)=o_o(i);
figure
subplot(3,1,1)
hold on
plot(t_s,o,'LineWidth',1.5);
plot(t_s, o_ref,'LineWidth',1.5, 'Color', 'g')
plot(t_s, o_o,'LineWidth',1.5, 'Color', 'y')
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
plot(t_s, ia_o, 'Color', 'b')
title('Salida del sistema i_a(t)')
xlabel('Tiempo [s]')
ylabel('Corriente [A]')
grid
hold off