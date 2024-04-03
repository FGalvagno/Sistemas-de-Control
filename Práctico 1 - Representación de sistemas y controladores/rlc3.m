% ITEM 3------------------------------------------------------------------------
[u,t] = gensig("square",0.08,0.16,0.00005);
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