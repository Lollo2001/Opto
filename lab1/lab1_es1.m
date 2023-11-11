clc
clear
close all

%% 1 - Valutazione n_eff(d)
sample = 10^4;  % numero di valori per la creazione dei vettori tramite la funzione linspace

n_cladding=1.44;
n_core=3.48;
c=3e8;
lambda=1.55*10^-6;
omega=2*pi*c/lambda;
mu=4*pi*10^-7;
m0=0;
m1=1;
m2=2;
k_0=2*pi/lambda;

n_eff=linspace(n_cladding,n_core-0.0001,sample);    % vettore di punti tra n_cladding e n_core (sample più grande di 1000 per avere una buona risuluzione)
gamma= k_0*(sqrt(n_eff.^2 -n_cladding^2));

k_x= k_0*sqrt(n_core^2-n_eff.^2);

d0=(2*atan(gamma./k_x)+m0*pi)./k_x;     % spessore del core per il modo 0
d1=(2*atan(gamma./k_x)+m1*pi)./k_x;     % spessore del core per il modo 1
d2=(2*atan(gamma./k_x)+m2*pi)./k_x;     % spessore del core per il modo 2

figure (1)
plot( d0*10^6, n_eff, d1*10^6, n_eff, d2*10^6, n_eff, 'LineWidth', 1.5)
xlabel('d [µm]');
ylabel('n_{eff}');
legend('TE_0','TE_1','TE_2');
title('n_{eff} in funzione di d');
xlim([0, 7]);   
ylim([n_cladding, n_core]);
grid on
%% punto 2
d=4.5e-7; % valore scelto per lo spessore del core per avere solo modo TE0 e TE1
[~, index_0] = min(abs(d0-d));  % trovo l'indice del dato del vettore d0 con il valore più vicino a d
[~, index_1] = min(abs(d1-d));  % trovo l'indice del dato del vettore d1 con il valore più vicino a d
n_eff_0 = n_eff(index_0);      % valore efficace per il modo 0
n_eff_1 = n_eff(index_1);      % valore efficace per il modo 1
beta_0 = k_0*n_eff_0;       % beta efficace per il modo 0
beta_1 = k_0*n_eff_1;       % beta efficace per il modo 1
k_x_0 = k_0*(sqrt(n_core^2-n_eff_0^2)); % k_x per il modo 0
k_x_1 = k_0*(sqrt(n_core^2-n_eff_1^2)); % k_x per il modo 1
gamma_0 = k_0*(sqrt(n_eff_0^2 -n_cladding^2));  % gamma per il modo 0
gamma_1 = k_0*(sqrt(n_eff_1^2 -n_cladding^2));  % gamma per il modo 1

% Ey
syms x;
% espressione analitica a tratti per il campo elettrico lungo l'asse y per il modo 0
Ey_0 = piecewise( x>0, exp(-gamma_0.*x), -d<=x<=0, (cos(k_x_0.*x) - (gamma_0/k_x_0)*sin(k_x_0.*x)), x<-d, (cos(k_x_0.*(-d)) - (gamma_0/k_x_0)*sin(k_x_0.*(-d)))*exp(gamma_0*(x+d)) );
figure(2)
hold on
grid on
fplot(Ey_0, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("E_y(x) di TE_0")
xlabel('x [m]')
ylabel('E_y [V/m]')
xlim([-2*d, d])

% espressione analitica a tratti per il campo elettrico lungo l'asse y per il modo 1
Ey_1 = piecewise( x>0, exp(-gamma_1.*x), -d<=x<=0, (cos(k_x_1.*x) - (gamma_1/k_x_1)*sin(k_x_1.*x)), x<-d, (cos(k_x_1.*(-d)) - (gamma_1/k_x_1)*sin(k_x_1.*(-d)))*exp(gamma_1*(x+d)) );
figure(3)
hold on
grid on
fplot(Ey_1, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("E_y(x) di TE_1")
xlabel('x [m]')
ylabel('E_y [V/m]')
xlim([-2*d, d])

% Hz
% espressione analitica a tratti per il campo magnetico lungo l'asse z per il modo 0
Hz_0 = diff(Ey_0, x)/(omega*mu);
figure(4)
figure(4)
hold on
grid on
fplot(Hz_0, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_z(x) di TE_0")
xlabel('x [m]')
ylabel('H_z [A/m]')
xlim([-2*d, d])

% espressione analitica a tratti per il campo magnetico lungo l'asse z per il modo 1
Hz_1 = diff(Ey_1, x)/(omega*mu);
figure(5)
hold on
grid on
fplot(Hz_1, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_z(x) di TE_1")
xlabel('x [m]')
ylabel('H_z [A/m]')
xlim([-2*d, d])

% Hx 
% espressione analitica a tratti per il campo magnetico lungo l'asse x per il modo 0
Hx_0 = -(beta_0*Ey_0)/(omega*mu);
figure(6)
hold on
grid on
fplot(Hx_0, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_x(x) di TE_0")
xlabel('x [m]')
ylabel('H_x [A/m]')
xlim([-2*d, d])

% espressione analitica a tratti per il campo magnetico lungo l'asse x per il modo 1
Hx_1 = -(beta_1*Ey_1)/(omega*mu);
figure(7)
hold on
grid on
fplot(Hx_1, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_x(x) di TE_1")
xlabel('x [m]')
ylabel('H_x [A/m]')
xlim([-2*d, d])
%}


%% punto 3
Gamma = ones(1,length(n_eff));  % inizializzo vettore in cui poi verrà salvato il valore di Gamma per ogni spessore del core

for i = 1:length(d0)     % per ogni valore del vettore d0 vado a calcolare il fattore Gamma e lo salvo nell'ononimo vettore
    syms x;
    % espressione analitica a tratti del campo elettrico lungo l'asse y per il modo 0
    Ey = piecewise( x>0, exp(-gamma(i).*x), -d0(i)<=x<=0, (cos(k_x(i).*x) - (gamma(i)/k_x(i))*sin(k_x(i).*x)), x<-d0(i), (cos(k_x(i).*(-d0(i))) - (gamma(i)/k_x(i))*sin(k_x(i).*(-d0(i))))*exp(gamma(i)*(x+d0(i))) );
    Gamma(i) = int(Ey^2, x, -d0(i), 0)/int(Ey^2, x, -inf, inf);    
end 

figure(8)
plot(d0*1e6, Gamma)
title("Γ_{TE0} in funzione dello spesso del core");
ylabel("Γ_{TE0}");
xlabel('d [µm]')
xlim([0, 1])
grid on



%% punto 4
alpha_core = 0.3;       % cm^-1
alpha_cladding = 0.15;  % cm^-1

alpha_modal = Gamma.*alpha_core + (1-Gamma).*alpha_cladding;    % alpha_modal in cm^-1

figure(9)
subplot(2, 1, 1)
plot(d0*1e6, alpha_modal)
title("Perdita modale in funzione dello spessore del core");
ylabel("\alpha_{modal} [{cm^-1}]");
xlabel('d [µm]')
xlim([0, 1])
grid on

alpha_modal_db = alpha_modal.*10*log10(exp(1));                 % alpha_modal_db in dB/cm
subplot(2, 1, 2)
plot(d0*1e6, alpha_modal_db) 
title("Perdita modale in funzione dello spessore del core");
ylabel("\alpha_{modal} [dB/cm]");
xlabel('d [µm]')
xlim([0, 1])
grid on


%% punto 5

len = 0.5;          % guida lunga 0.5 cm
insertion_loss = alpha_modal_db*len; 
figure(10)
plot(d0*1e6, insertion_loss)
title("Insertion loss in funzione dello spessore del core");
ylabel("Insertion loss [dB]");
xlabel('d [µm]')
xlim([0, 1])
grid on

%% punto 6
[~, index_Gamma] = min(abs(Gamma-0.75));  % trovo l'indice del dato del vettore Gamma con il valore più vicino a 0.75
thickness = d0(index_Gamma);    % ricavo lo spessore del core per il modo TE0 e con gamma pari a 0.75

%% punto 7
lambda_min = 850e-9;    % valore minimo per la lunghezza d'onda da raffigurare sul grafico
lambda_max = 1620e-9;   % valore massimo per la lunghezza d'onda da raggigurare sul grafico
w_0 = 2.*atan(sqrt(n_eff.^2 - n_cladding^2)./sqrt(n_core^2 - n_eff.^2)).*c./(thickness.*sqrt(n_core^2 - n_eff.^2));         % pulsazione per il modo 0
w_1 = (2.*atan(sqrt(n_eff.^2 - n_cladding^2)./sqrt(n_core^2 - n_eff.^2))+pi).*c./(thickness.*sqrt(n_core^2 - n_eff.^2));    % pulsazione per il modo 1

subplot(1,2,1)
plot(w_0, n_eff, w_1, n_eff)
title("n_{eff} in funzione di omega");
ylabel("n_{eff}");
xlabel('omega [rad/s]')
xlim([2*pi*c/lambda_max, 2*pi*c/lambda_min])    % limito il grafico lungo l'asse x in base ai valori di lambda_min e lambda_max dati
ylim([n_cladding-0.1, n_core+0.1])              % limito il grafico lungo l'asse y per una migliore visualizzazione
yline(n_core, 'k--')        % inserisco sul grafico una retta orizzonatale che raffigura il valore di n_core
yline(n_cladding, 'k-.')    % inserisco sul grafico una retta orizzonatale che raffigura il valore di n_cladding
legend('TE_0', 'TE_1', 'n_{core}', 'n_{cladding}')
grid on

lambda_0=2*pi*c./w_0;   % dal valore di pulsazione per il modo 0 calcolato prima, calcolo la lunghezza d'onda per il modo 0
lambda_1=2*pi*c./w_1;   % dal valore di pulsazione per il modo 1 calcolato prima, calcolo la lunghezza d'onda per il modo 1
subplot(1,2,2)
plot(lambda_0*1e9, n_eff, lambda_1*1e9, n_eff)
title("n_{eff} per TE_0 in funzione di lambda");
ylabel("n_{eff}");
xlabel('lambda [nm]')
xlim([lambda_min*1e9, lambda_max*1e9])  % limito il grafico lungo l'asse x in base ai valori di lambda_min e lambda_max dati
ylim([n_cladding-0.1, n_core+0.1])      % limito il grafico lungo l'asse y per una migliore visualizzazione
yline(n_core, 'k--')        % inserisco sul grafico una retta orizzonatale che raffigura il valore di n_core
yline(n_cladding, 'k-.')    % inserisco sul grafico una retta orizzonatale che raffigura il valore di n_cladding
legend('TE_0', 'TE_1', 'n_{core}', 'n_{cladding}')
grid on




















