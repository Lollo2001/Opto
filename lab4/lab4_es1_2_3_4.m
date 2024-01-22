clc
clear 
close all

%% es1

lambda=1550e-9; % metri
n_eff=3.86; 
L= 5e-3; % metri
V = 0:0.01:5; % V
d_phi=1.732*sqrt(V);    % variazione di fase in funzione della tensione applicata
% figure()
% hold on
% plot(V,d_phi,'LineWidth',1.5)
% grid on
T_loss=(-1.5*V+21.5);   % Trasmission loss [dB] in funzione della tensione applicata
alfa=T_loss/(L*100*4.34);   % equazione per calcolare le perdite (alfa) [cm^-1] a partire dal coefficiente di trasmissione
figure()
plot(V,alfa,'LineWidth',1.5)    % grafico perdite dovute all'attenuazione della guida in funzione della tensione applicata
title("Perdite della guida in funzione della tensione inversa applicata alla giunzione p-i-n")
xlabel("Tensione [V]")
ylabel("α_{TOT} [cm^{-1}]")
grid on
%phi=2*pi*L*n_eff/lambda + d_phi; 
d_neff=d_phi*lambda/(2*pi*L);   % data la variazione di fase calcolo la variazione di indice efficace
%neff_tot=n_eff+d_neff;
figure()
plot(V,d_neff,'LineWidth',1.5)  % grafico dell'indice di rifrazione in funzione della tensione applicata
title("Indice di rifrazione in funzione della tensione inversa applicata alla giunzione p-i-n")
xlabel("Tensione [V]")
ylabel("Δn_{eff}")
grid on

%% es2
lambda=linspace(1545,1565,1e4)*1e-9;    % metri
d_L = 1e-4;     % Il ramo del modulatore in cui NON applichiamo tensione è lungo L+dL=5.1mm
alfa_0 = 21.5/(L*4.34); % in metri
V=0:5;
figure()
grid on
for i = 1:6    % simulazione del coefficente di Trasmissione al variare della lunghezza d'onda a varie tensione [da 0 a 5V]
    d_phi=1.732*sqrt(V(i));     % dal grafico
    d_neff=d_phi.*lambda./(2*pi*L);   % variazione dell'indicide di rifrazione efficace in funzione della tensione applicata
    T_loss=(-1.5*V(i)+21.5);   % dal grafico
    alfa=T_loss/(L*4.34);       % in metri
    
    campo1 = ( exp(-1i.*2.*pi.*L.*(n_eff+d_neff)./lambda)*exp(-alfa*L/2) )/ sqrt(2); 
    campo2 = (exp(-1i.*2.*pi.*(L+d_L).*n_eff./lambda)*exp(-alfa_0*(L+d_L)/2) )/ sqrt(2);

    T = abs((campo1 + campo2)./sqrt(2)).^2;
    T_db=10*log10(T); 
    hold on 
    plot(lambda*1e9,T_db,'LineWidth',1.5) 
end
title("Spettro di trasmissione in funzione di lambda a diversi valori di tensione applicata")
xlabel("Lambda [nm]")
ylabel("Spettro di trasmissione [dB]")
legend("V=0", "V=1", "V=2", "V=3", "V=4", "V=5")

%% es3

% Cerco la lunghezza d'onda (lambda_in) corrispondente ad un massimo di T per V=0
lambda=linspace(1545,1565,1e4)*1e-9;    % metri
d_L = 1e-4;     % Il ramo del modulatore in cui NON applichiamo tensione è lungo L+dL=5.1mm
alfa_0 = 21.5/(L*4.34); % in metri
campo1 = ( exp(-1i.*2.*pi.*L.*(n_eff)./lambda)*exp(-alfa_0*L/2) )/ sqrt(2); 
campo2 = (exp(-1i.*2.*pi.*(L+d_L).*n_eff./lambda)*exp(-alfa_0*(L+d_L)/2) )/ sqrt(2);
T = abs((campo1 + campo2)./sqrt(2)).^2;
[~, indice] = max(T);    % cerco il valore dell'indice del primo massimo di T
lambda_in = lambda(indice)  % associo il corrispettivo valore di lambda

% Coefficiente di trasmissione in funzione della tensione per lambda=lambda_in
V=linspace(0,6,1e6);    % vettore delle tensioni
d_phi=1.732.*sqrt(V);     % dal grafico
d_neff=d_phi.*lambda_in./(2*pi*L);   % variazione dell'indicide di rifrazione efficace in funzione della tensione applicata
T_loss=(-1.5*V+21.5);   % dal grafico
alfa=T_loss/(L*4.34);       % in metri
alfa_0 = alfa(1);       % alfa quando V=0
campo1 = ( exp(-1i.*2.*pi.*L.*(n_eff+d_neff)./lambda_in).*exp(-alfa*L/2) )/ sqrt(2); 
campo2 = (exp(-1i.*2.*pi.*(L+d_L).*n_eff./lambda_in).*exp(-alfa_0*(L+d_L)/2) )/ sqrt(2);
T = abs((campo1 + campo2)./sqrt(2)).^2;     % coefficiente di trasmissione in lineare
T_db=10*log10(T);       % coefficiente di trasmissione in dB
figure()
plot(V,T_db,'LineWidth',1.5) 
grid on
title("Coefficiente di trasmissione in funzione della tensione applicata")
xlabel("Tensione [V]")
ylabel("Coefficiente di trasmissione [dB]")

% Cerco V_pigreco (minimo del grafico del coefficiente di trasmissione appena fatto)
[~, indice] = min(T_db);
Vpi = V(indice)     % V_pigreco

% Cerco Vbias che massimizza ER, ovvero dove la derivata di T_db è maggiore
derivata = diff(T_db);
[~, indice] = min(derivata);    % voglio trovare l'indice del punto che ha la derivata con segno negativo e valore più grande possibile
Vbias = V(indice+1)     % aggiungo 1 perché la derivata mi restituisce N-1 valori
Vpp = 1;    % V
V0 = Vbias+Vpp/2;
V1 = Vbias-Vpp/2;
[~, indiceV0] = min(abs(V0-V));
T_V0 = T(indiceV0);
[~, indiceV1] = min(abs(V1-V));
T_V1 = T(indiceV1);
ER = 10*log10(T_V1/T_V0)

% Calcolo il Transmission Loss in funzione della tensione di bias al livello 1
Vbias=linspace(0+Vpp/2,6+Vpp/2,1e6);    % vettore delle tensioni Vbias
T_loss=(-1.5*(Vbias-Vpp/2)+21.5);   % dal grafico ma bisogna modificare la tensione per avere le perdite del livello 1 in funzione di Vbias
figure()
plot(Vbias,T_loss,'LineWidth',1.5) 
grid on
title("Perdite di trasmissione del livello 1 in funzione della tensione applicata Vbias")
xlabel("Tensione Vbias [V]")
ylabel("Perdite di trasmissione [dB]")


%% es4
f = linspace(0.1,1e12,100000);
RL = 37;
RG = 50;
Z0 = 37;
ZL = RL;
L = 5e-3;    % m
omega = 2*pi.*f;
n_o = 3.75;
n_m = 3.86;
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
figure()
semilogx(f, m_dB, 'LineWidth',1.5)
yline(m_dB(1) - 3);
grid on
title("m(f) per una guida senza perdite e adatatta")
xlabel("frequenza [Hz]")
ylabel("m(f) [dB]")
xlim([0,1e12])
%legend("m(f)", "m(0) - 3dB")

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
% m= abs(Fpiu./exp(gamma.*L));  % espressione di m semplificato
m_dB = 10.*log(m);
hold on 
semilogx(f, m_dB, 'LineWidth',1.5)
%yline(m_dB(1) - 3);
title("m(f) per una guida senza perdite e adatatta")
xlabel("frequenza [Hz]")
ylabel("m(f) [dB]")
legend("m(f) senza perdite", "m(0) - 3dB", "m(f) con perdite")

[~, indice_B] = min(abs(m_dB - (m_dB(1) - 3)));
f_meno3dB_B = f(indice_B)   % frequenza a -3dB

%% punto c
%energia trasmissione
E = (0.5^2)/(RG+RL); %energia al secondo
Ebit = E/f_meno3dB_B    % energia per bit
