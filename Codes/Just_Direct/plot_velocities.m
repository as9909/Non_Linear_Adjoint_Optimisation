close all;clear;clc
%%
U=importdata('U_vel_stretched_grid_steady_stratified_viscosity.txt');
% U2=importdata('U2_vel.txt');

yf=importdata('yf_grid_stretched.txt');
% 
figure,plot(U(1,2:end-1),yf(2:end-1),'k-','LineWidth',3)
hold on
plot(U(5,2:end-1),yf(2:end-1),'-','LineWidth',2)
plot(U(10,2:end-1),yf(2:end-1),'-','LineWidth',2)
plot(U(20,2:end-1),yf(2:end-1),'-','LineWidth',2)
plot(U(50,2:end-1),yf(2:end-1),'-','LineWidth',2)
plot(U(200,2:end-1),yf(2:end-1),'-','LineWidth',2)
plot(U(500,2:end-1),yf(2:end-1),'-','LineWidth',2)
plot(U(1000,2:end-1),yf(2:end-1),'-','LineWidth',2)
 plot(U(3000,2:end-1),yf(2:end-1),'-','LineWidth',2)
% plot(U(4500,2:end-1),yf(2:end-1),'-','LineWidth',2)
xlabel('$u$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
set(gca,'ticklabelInterpreter','latex','FontSize',20)
grid on
K1=0.6;
u_theory=-2/K1*(1+coth(K1)+(yf-coth(K1)).*exp(K1*(1+yf)));
%u_theory=3.0/2*1*(1-(2*yf/2).^2);
plot(u_theory(2:end-1),yf(2:end-1),'ko','LineWidth',2)
set(gca,'LooseInset',get(gca,'TightInset'))
print('stretched_grid_steady_stratified_viscosity.eps','-depsc')
for ii=1:length(U(:,1))
err(ii)=trapz(yf(2:end-1),abs(u_theory(2:end-1)-U(ii,(2:end-1))));
end
figure,semilogy([0:length(err)-1],err,'k-','LineWidth',2)
xlabel('Time','Interpreter','latex')
ylabel('Error','Interpreter','latex')
grid on
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'ticklabelInterpreter','latex','FontSize',20)
print('Error_stretched_grid_steady_stratified_viscosity.eps','-depsc')
%figure,plot(T(1500,2:end-1),'--','LineWidth',2)
%{
AA=U-u_theory;
A_c=max(abs(AA(:,2:end-1)),[],2);
figure,semilogy(A_c)
% figure,plot((U(round(end/2),2:end-1)-U(end,2:end-1)))
% 
% figure,plot(yf,mu(1,:),'.')
% figure,plot(U2(1,2:end-1))
% hold on
% plot(U2(100,2:end-1),'--')
% plot(U2(300,2:end-1),'-.')
% title('u')
% 
% figure,plot(pr1(1,2:end-1))
% hold on
% plot(pr1(100,2:end-1))
% plot(pr1(300,2:end-1))
% title('pr')

% figure,plot(chk_Poss_RHS(1,2:end-1))
% hold on
% plot(chk_Poss_RHS(100,2:end-1))
% plot(chk_Poss_RHS(300,2:end-1))
% title('check pressure')


% 
% figure,plot(U(1,2:end-1)- U2(1,2:end-1))
% hold on
% plot(U(100,2:end-1)- U2(100,2:end-1))
% plot(U(300,2:end-1)- U2(300,2:end-1))
% title('u at diff point')
% 
% figure,plot(y(2:end-1),V(1,2:end-1))
% hold on
% plot(y(2:end-1),V(10,2:end-1))
% plot(y(2:end-1),V(330,2:end-1))
% title('v')
% 
% 
% figure,plot(W(1,2:end-1))
% hold on
% plot(W(10,2:end-1))
% plot(W(330,2:end-1))
% title('w')
% 
% 
% figure,plot(mu(1,2:end-1))
% hold on
% plot(mu(10,2:end-1))
% plot(mu(330,2:end-1))
% title('mu')
%}