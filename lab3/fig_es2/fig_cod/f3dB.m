clear
close all
clc

L_a=9e-2; % cm
r1=sqrt(0.3);
r2=sqrt(0.9);
alfa_i=6; % cm-1
alfa_m=log(1/(r1*r2))/L_a; % cm-1
c=3e10; % cm/s
eta_g=4.2;
v_g=c/eta_g; % cm/s
gain=1e-16; % cm3
a0=1.5e-16; % cm2
I_th=10; % mA circa

tau_p=1/((alfa_i+alfa_m)*v_g); % s
I_bias=linspace(I_th,10*I_th,1e5); % mA
f_max=sqrt(2)/(2*pi*tau_p*1e9); % GHz
f_max_gain=f_max*(1/(1+(alfa_m+alfa_i)*gain/a0)); % GHz

I=[1.05 2 6 10]*I_th; % mA
f=[0.69027 2.96275 6.43918 8.52135]; % GHz
f_gain=[0.68438 2.14182 2.17815 2.09905]; % GHz

figure(1)
plot(I,f,'ro-','LineWidth',1.5)
hold on
yline(f_max,'b-','LineWidth',1.5)
xlabel("I_{bias} [mA]",FontSize=14)
ylabel("f_{-3dB} [GHz]",FontSize=14)
legend("f_{-3dB}","f_{-3dB,max}",FontSize=14)

figure(2)
plot(I,f_gain,'ro-','LineWidth',1.5)
hold on
yline(f_max_gain,'b-','LineWidth',1.5)
xlabel("I_{bias} [mA]",FontSize=14)
ylabel("f_{-3dB} con ε [GHz]",FontSize=14)
legend("f_{-3dB}","f_{-3dB,max} con ε",FontSize=14)