clc
clear 
close all

%% es1

lambda=1550e-7; 
n_eff=3.86; 
L=0.5; 
V = 0:0.01:5; 
d_phi=1.732*sqrt(V); 
% figure(1)
% hold on
% plot(V,d_phi,'LineWidth',1.5)
% grid on
T_loss=(-1.5*V+21.5); 
alfa=T_loss/(L*4.34); 
figure(2)
plot(V,alfa,'LineWidth',1.5) 
grid on
phi=2*pi*L*n_eff/lambda + d_phi; 
d_neff=d_phi*lambda/(2*pi*L); 
neff_tot=n_eff+d_neff;
figure(3)
plot(V,d_neff,'LineWidth',1.5) 
grid on

%% es2

lambda=linspace(1545,1565,1e4)*1e-7; 
d_L=100e-4; 
V=0:5;
alfa_0=21.5/(L*4.34); 
for i = 1:6
    d_phi=1.732*sqrt(V(i)); 
    d_alfa=-1.5*V(i)/(4.34*L); 
    d_neff=d_phi.*lambda./(2*pi*L); 
    T=exp(-alfa_0*L)/4*abs(exp(-1i*2*pi*d_neff*L./lambda)*exp(-d_alfa*L/2)+exp(-alfa_0*d_L/2)*exp(-1i*2*pi*n_eff*d_L./lambda)).^2; 
    T_db=10*log10(T); 
    figure(4)
    hold on 
    plot(lambda*1e7,T_db,'LineWidth',1.5) 
end

%% es3

lambda=1550e-7;
V=linspace(0,6,1e6); 
alfa_0=21.5/(L*4.34);
d_phi=1.732.*sqrt(V); 
d_alfa=-1.5.*V./(4.34*L); 
d_neff=d_phi.*lambda/(2*pi*L); 
T=exp(-alfa_0.*L)./4.*abs(exp(-1i*2*pi.*d_neff.*L/lambda).*exp(-d_alfa.*L/2)+exp(-alfa_0.*d_L/2).*exp(-1i*2*pi.*n_eff.*d_L/lambda)).^2; 
T_db=10*log10(T); 
figure(5)
plot(V,T_db,'LineWidth',1.5) 
[min_T,index_T]=min(T_db); 
V_pi=V(index_T); 
V_pp=1; 
V_p=V_pp/2; 
