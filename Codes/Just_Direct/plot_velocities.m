close all;clear;clc
%%
U=importdata('U_vel.txt');
U2=importdata('U2_vel.txt');
V=importdata('V_vel.txt');
W=importdata('W_vel.txt');
mu=importdata('mu.txt');
pr1=importdata('pressure1.txt');
pr2=importdata('pressure2.txt');
pr3=importdata('pressure3.txt');
% pr4=importdata('pressure4.txt');
% chk_Poss_RHS=importdata('chk_Poss_RHS.txt');

y=importdata('y_grid.txt');
yf=importdata('yf_grid.txt');
% 
figure,plot(yf(2:end-1),U(1,2:end-1))
hold on
plot(yf(2:end-1),U(5,2:end-1),'-','LineWidth',2)
plot(yf(2:end-1),U(10,2:end-1),'-','LineWidth',2)
plot(yf(2:end-1),U(20,2:end-1),'-','LineWidth',2)
plot(yf(2:end-1),U(50,2:end-1),'-','LineWidth',2)
plot(yf(2:end-1),U(500,2:end-1),'-','LineWidth',2)
plot(yf(2:end-1),U(1500,2:end-1),'--','LineWidth',2)
title('u')
%{
u_theory=3.0/2*1*(1.0-(2.0*yf).^2);
AA=U-u_theory;
A_c=max(abs(AA(:,2:end-1)),[],2);
figure,semilogy(A_c)
%}
% figure,plot(yf(2:end-1),(U(round(end/2),2:end-1)-U(end,2:end-1)))
% 
% figure,plot(yf,mu(1,:),'.')
% figure,plot(yf(2:end-1),U2(1,2:end-1))
% hold on
% plot(yf(2:end-1),U2(100,2:end-1),'--')
% plot(yf(2:end-1),U2(300,2:end-1),'-.')
% title('u')
% 
% figure,plot(yf(2:end-1),pr1(1,2:end-1))
% hold on
% plot(yf(2:end-1),pr1(100,2:end-1))
% plot(yf(2:end-1),pr1(300,2:end-1))
% title('pr')

% figure,plot(yf(2:end-1),chk_Poss_RHS(1,2:end-1))
% hold on
% plot(yf(2:end-1),chk_Poss_RHS(100,2:end-1))
% plot(yf(2:end-1),chk_Poss_RHS(300,2:end-1))
% title('check pressure')


% 
% figure,plot(yf(2:end-1),U(1,2:end-1)- U2(1,2:end-1))
% hold on
% plot(yf(2:end-1),U(100,2:end-1)- U2(100,2:end-1))
% plot(yf(2:end-1),U(300,2:end-1)- U2(300,2:end-1))
% title('u at diff point')
% 
% figure,plot(y(2:end-1),V(1,2:end-1))
% hold on
% plot(y(2:end-1),V(10,2:end-1))
% plot(y(2:end-1),V(330,2:end-1))
% title('v')
% 
% 
% figure,plot(yf(2:end-1),W(1,2:end-1))
% hold on
% plot(yf(2:end-1),W(10,2:end-1))
% plot(yf(2:end-1),W(330,2:end-1))
% title('w')
% 
% 
% figure,plot(yf(2:end-1),mu(1,2:end-1))
% hold on
% plot(yf(2:end-1),mu(10,2:end-1))
% plot(yf(2:end-1),mu(330,2:end-1))
% title('mu')

