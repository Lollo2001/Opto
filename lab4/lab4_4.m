clc
clear
close all

f = linspace(0.1,1e12,100000);
RL = 37;
RG = 50;
Z0 = 37;
ZL = RL;
L = 5e-3;    % m
omega = 2*pi.*f;
n_m = 3.75;
n_o = 3.86;
c = 3e8;    % m/s
alpha_dB = 10.8;    %dB/cm

%% punto a
alpha = 0;   %attenuazione al metro
gamma = alpha + 1i.*omega.*n_m./c;
Zin = Z0.*(ZL+Z0.*atan(gamma.*L))./(Z0+ZL.*atan(gamma.*L));
u_piu = alpha.*L + 1i.*omega./c.*(n_m-n_o).*L; 
u_meno = -alpha.*L + 1i.*omega./c.*(-n_m-n_o).*L; 
Fpiu = (1-exp(u_piu))./u_piu;
Fmeno = (1-exp(u_meno))./u_meno;

m = ((RL+RG)./RL).*abs(Zin./(Zin+RG)).*abs(((ZL+Z0).*Fpiu + (ZL-Z0).*Fmeno)./((ZL+Z0).*exp(gamma.*L) + (ZL-Z0).*exp(-gamma.*L)));
m_dB = 10.*log(m);
semilogx(f, m_dB)
yline(m_dB(1) - 3);
grid on
title("m(f) per una guida senza perdite e adatatta")
xlabel("frequenza [Hz]")
ylabel("m(f) [dB]")
xlim([0,1e12])

[~, indice_A] = min(abs(m_dB - (m_dB(1) - 3)));
f_meno3dB_A = f(indice_A)   % frequenza a -3dB

%% punto b
% alpha è espresso come attenuazione di potenza o di tensione?
alpha_10ghz = 10^(alpha_dB/20) *100; %attenuazione al metro a 10GHz
alpha0 = alpha_10ghz/(1e10)^0.5;    %attenuzione a 0Hz
alpha = alpha0.*f.^0.5; %attenuazione in funzione della frequenza
gamma = alpha + 1i.*omega.*n_m./c;
Zin = Z0.*(ZL+Z0.*atan(gamma.*L))./(Z0+ZL.*atan(gamma.*L));
u_piu = alpha.*L + 1i.*omega./c.*(n_m-n_o).*L; 
u_meno = -alpha.*L + 1i.*omega./c.*(-n_m-n_o).*L; 
Fpiu = (1-exp(u_piu))./u_piu;
Fmeno = (1-exp(u_meno))./u_meno;

m = ((RL+RG)./RL).*abs(Zin./(Zin+RG)).*abs(((ZL+Z0).*Fpiu + (ZL-Z0).*Fmeno)./((ZL+Z0).*exp(gamma.*L) + (ZL-Z0).*exp(-gamma.*L)));
m_dB = 10.*log(m);
hold on  
semilogx(f, m_dB)
yline(m_dB(1) - 3);
title("m(f) per una guida senza perdite e adatatta")
xlabel("frequenza [Hz]")
ylabel("m(f) [dB]")

[~, indice_B] = min(abs(m_dB - (m_dB(1) - 3)));
f_meno3dB_B = f(indice_B)   % frequenza a -3dB