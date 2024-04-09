% ITEM 1------------------------------------------------------------------------
clear all; close all; clc
%Generamos señal de entrada:
%Onda cuadrada alterna de 2ms de periodo y 12V de amplitud


R=47; L=1e-06; C=100e-09; %Parámetros del circuito RLC
A1 = [-R/L -1/L; 1/C 0];
B1 = [1/L; 0];
C1 = [R 0];
C2 = [0 1]; 
D1 = [0];

%Obtenemos la FdT
[N,D] = ss2tf(A1,B1,C1,D1);
FdT = tf(N,D)

%Estudiamos la dinámica en base a los autovalores de A
eig_val = eig(A1)

%Tiempo de itegración y simulación
tr=log(0.95)/eig_val(1);

%Tomamos ti = 1e-10 (10 veces mas pequeño)
ti = 1e-10

%Tomamos ts = 5e-03 debido a que buscamos simular la respuesta con una señal
%que conmuta
ts = 5e-03 

N = ts/ti;

%Señal de entrada y vector de tiempo
t = linspace (0, ts, N);
u = linspace (0, 0, N);

%Funciones i_c, v_c 

i_c=zeros(1,N);
v_c=zeros(1,N);


x = [0 0]'; %Cond. iniciales nulas
vin = 12;
ij = 0;

for ii=1:N-1
    ij=ij+ti;
    if (ij>=1e-3)
        ij=0;
        vin=vin*-1;
    end
    u(ii)=vin;
    xp=A1*x+B1*u(ii);
    x=x+xp*ti;
    Y=C1*x;
    y(ii+1)=Y(1);
    i_c(ii+1)=x(1);
    v_c(ii+1)=x(2);
end

figure('Name','1.1')
 subplot(3,1,1);plot(t,i_c);grid on; title('Corriente, I'); 
 subplot(3,1,2);plot(t,v_c);grid on; title('Tensión del capacitor, V_c');
 subplot(3,1,3);plot(t,u);grid on; title('Tensión de entrada, V_e');
 