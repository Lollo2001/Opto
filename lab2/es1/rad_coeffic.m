clear
close all
clc

A=3.57e8;
B=0.8e-10;
C=3.5e-30;
N=linspace(0,1e20,1e6);

rad_cf=(B.*N)./(A+B.*N+C.*N.^2);
figure(1)
plot(N,rad_cf,'LineWidth',2)
title('Radiactive coefficiency along Carrier Density')
xlabel('N [cm^{-3}]')
ylabel('η_r')
grid on

max_values_A=zeros([1 4]);
max_values_B=zeros([1 4]);
max_values_C=zeros([1 4]);

for k = 1:4
    rad_cf_A=(B.*N)./(A*k+B.*N+C.*N.^2);
    rad_cf_B=(k*B.*N)./(A+k*B.*N+C.*N.^2);
    rad_cf_C=(B.*N)./(A+B.*N+k*C.*N.^2);
    max_values_A(k)=max(rad_cf_A);
    max_values_B(k)=max(rad_cf_B);
    max_values_C(k)=max(rad_cf_C);
    % figure(2)
    % hold on
    % plot(N,rad_cf_A,'LineWidth',2)
    % figure(3)
    % hold on
    % plot(N,rad_cf_B,'LineWidth',2)
    % figure(4)
    % hold on
    % plot(N,rad_cf_C,'LineWidth',2)
end