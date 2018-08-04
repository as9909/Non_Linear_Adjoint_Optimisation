

cd '/home/ritabrata/Desktop/Non_Linear_Adjoint_Optimisation-master/Codes/Just_Direct'
cd '/home/ritabrata.thakur/Desktop/Just_Direct'
u_64points_075stretch = importdata('U_vel.txt');
yf_64points_075stretch  = importdata('yf_grid.txt');



figure
plot(yf_64points_075stretch(2:end-1), u_64points_075stretch(1,2:end-1), '*')
hold on
plot(yf_64points_075stretch(2:end-1), u_64points_075stretch(1000,2:end-1))
legend('theroretical', 'After 1000 time steps')
title('viscosity = (1 - 0.25*T), medium stretching, more NY')


for i = 2:length(U)
error_64Ny(i) = trapz(yf_64points_075stretch(2:end-1), abs(u_64points_075stretch(i,2:end-1)-u_64points_075stretch(1,2:end-1)));
%max(u(1,2:end-1) - u())
end

clear error_new
for i = 2:length(U)
error_new(i) = max(abs(U(i,2:end-1)-U(1,2:end-1)));
%max(u(1,2:end-1) - u())
end

figure 
plot(error_new)
title('0.75 stretching, NY = 128')

figure 
plot(error_64Ny_2)
xlabel('time steps')
ylabel('Error')
title('maximum error at each time step, 0.75 stretching, NY = 128')

