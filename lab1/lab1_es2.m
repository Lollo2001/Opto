close all
clear
clc

%% punto 1 (impulsi di ingresso e spettri relativi con tau_0 variabili)

tau_0=[1e-9 5e-10 1e-10 1e-11]; % costanti di tempo [s]
N=1e3; % fattore di campionamento
M=70; % fattore per l'ampiezza del segnale
c=3e8; % velocità della luce [m/s]
lambda_w=1625e-9; % lunghezza d'onda in esame
f_0=c/lambda_w; % frequenza della portante

for j = 1:length(tau_0) % ciclo per gestire ogni costante di tempo
    t_end=M*tau_0(j); % limite per il vettore tempi
    dt=tau_0(j)/N; % tempo di campionamento
    t=(-t_end:dt:t_end); % vettore tempi [s]
    lung=length(t); % lunghezza vettore tempi
    f=linspace(-0.5/dt,0.5/dt,lung)+f_0; % vettore frequenze [Hz]
    p_t=exp(-(t.^2/(2*tau_0(j)^2))); % segnale di ingresso di tipo gaussiano 
    P_f=fftshift(fft(p_t))/lung; % trasformata di Fourier per trovare lo spettro del segnale d'ingresso
    lambda=c./f; % lunghezze d'onda [m]
    figure(1)
    hold on  
    plot(t*1e9,p_t,'LineWidth',0.5)
    figure (2)
    hold on               %% plot ancora molto dubbio
    plot(lambda,abs(P_f),'LineWidth',0.5)
end

figure (1)
xlim([-5 5])
title('Segnale gaussiano di ingresso','FontSize',22)
xlabel('t [ns]','FontSize',16)           % plot dell'impulso d'ingresso
ylabel('|P_t|','FontSize',16)
legend('τ_0=1 ns','τ_0=0.5 ns','τ_0=0.1 ns','τ_0=10 ps','Location','northeast')
grid on

figure (2)
xlim([lambda_w-0.7e-9 lambda_w+0.7e-9])
ylim([0 max(abs(P_f))])
title('Spettro ottico del segnale gaussiano','FontSize',22)
xlabel('λ [m]','FontSize',16)           % plot dello spettro dell'ingresso
ylabel('|P_f|','FontSize',16)
legend('τ_0=1 ns','τ_0=0.5 ns','τ_0=0.1 ns','τ_0=10 ps','Location','southeast')
grid on

%% punto 2 (costante di propagazione β lungo ω)

D=22e-6; % valore massimo della dispersione della fibra
w_0=2*pi*f_0; % frequenza della portante
n_eff=1.463; % indice di rifrazione della fibra
n_g=1.4682; % indice di rifrazione di gruppo

w=linspace(0,10*w_0,1e4); % vettore pulsazioni [rad/s]

Beta_0=w_0/c*n_eff; % termine di ordine 0 dello sviluppo di taylor (valore alla pulsazione w_0)
Beta_0_primo=n_g/c; % termine di ordine 1 dello sviluppo di taylor (valore che determina il ritardo di gruppo)
Beta_0_secondo=-D*lambda_w^2/(2*pi*c); % termine di ordine 2 dello sviluppo di taylor (valore che determina la perdita della fibra)

Beta=Beta_0+Beta_0_primo.*(w-w_0)+0.5*Beta_0_secondo.*(w-w_0).^2; % costante di propagazione

figure(3)
plot(w,Beta,'LineWidth',0.5)
ylim([0 5e7])                     % plot della costante di propagazione
title('Costante di propagazione lungo ω','FontSize',22)
xlabel('ω [rad/s]','FontSize',16)
ylabel('β','FontSize',16)
grid on

%% punto 3 (impulsi e spettri di uscita con tau_0 e L variabili)

L=[1 10 1000 80000]; % lunghezze [m]
tau_0=[1e-9 1e-11]; % costanti di tempo [s]
alpha_db=0.00023; % attenuazione della fibra [dB/m]
N=1e3; % fattore di campionamento
M=70; % fattore per l'ampiezza del segnale

for k = 1:length(tau_0) % ciclo per gestire ogni costante di tempo
    for j = 1:length(L) % ciclo per gestire ogni lunghezza
        t_end=M*tau_0(k); % limite per il vettore tempi
        dt=tau_0(k)/N; % tempo di campionamento
        t=(-t_end:dt:t_end); % vettore tempi [s]
        lung=length(t); % lunghezza vettore tempi
        w=2*pi.*(linspace(-0.5/dt,0.5/dt,lung))+w_0; % vettore pulsazioni [rad/s]
        E_t=(tau_0(k)/sqrt(tau_0(k)^2+1i*Beta_0_secondo*L(j))).*exp(-(t.^2).*(tau_0(k)^2)./(2*(tau_0(k)^4+(Beta_0_secondo*L(j))^2))).*10^(-alpha_db/20*L(j)); % formula analitica dell'impulso d'uscita
        figure(4+2*(k-1))
        hold on
        subplot(1,2,1); plot(t*1e9,abs(E_t),'LineWidth',0.5)
        E_w=fftshift(fft(E_t))/lung; % trasformata di Fourier per trovare lo spettro del segnale d'uscita
        figure(5+2*(k-1))
        hold on
        subplot(1,2,1); plot(w,abs(E_w),'LineWidth',0.5)
    end
    p_t=exp(-(t.^2/(2*tau_0(k)^2))); % segnale di ingresso di tipo gaussiano 
    figure(4+2*(k-1))
    hold on     % plot dell'impulso d'ingresso
    subplot(1,2,1); plot(t*1e9,p_t,'LineWidth',0.5);
    P_w=fftshift(fft(p_t))/lung; % trasformata di Fourier per trovare lo spettro del segnale d'ingresso
    figure(5+2*(k-1))
    hold on      % plot dello spettro dell'ingresso
    subplot(1,2,1); plot(w,abs(P_w),'LineWidth',0.5)
end

for k = 1:length(tau_0) % ciclo per gestire ogni costante di tempo
    t_end=M*tau_0(k); % limite per il vettore tempi
    dt=tau_0(k)/N; % tempo di campionamento
    t=(-t_end:dt:t_end); % vettore tempi [s]
    lung=length(t); % lunghezza vettore tempi
    w=2*pi.*(linspace(-0.5/dt,0.5/dt,lung))+w_0; % vettore pulsazioni [rad/s]
    p_t=exp(-(t.^2/(2*tau_0(k)^2)));
    p_w=fftshift(fft(p_t))/lung; 
    for j = 1:length(L) % ciclo per gestire ogni lunghezza
        E_w=p_w.*exp(-0.5i.*Beta_0_secondo.*(w-w_0).^2.*L(j)).*10^(-alpha_db*L(j)/20); % spettro d'uscita
        E_t=ifft(E_w)*lung;  % impulso d'uscita  
        figure(4+2*(k-1))
        hold on      % plot dell'impulso d'uscita
        subplot(1,2,2); plot(t*1e9,abs(E_t),'LineWidth',0.5)
        figure(5+2*(k-1))
        hold on      % plot dello spettro d'uscita
        subplot(1,2,2); plot(w,abs(E_w),'LineWidth',0.5)
    end
    figure(4+2*(k-1))
    hold on     % plot dell'impulso d'ingresso
    subplot(1,2,2); plot(t*1e9,p_t,'LineWidth',0.5);
    figure(5+2*(k-1))
    hold on      % plot dello spettro d'ingresso
    subplot(1,2,2); plot(w,abs(p_w),'LineWidth',0.5)
end

figure(4)
subplot(1,2,1)
xlim([-5 5])
title('|E_t| con τ_0=1 ns analitica','FontSize',22)
xlabel('t [ns]','FontSize',16)      % plot dell'impulso d'ingresso analitico con tau_0= 1 ns
ylabel('|E_t|','FontSize',16)
legend('L=1 m','L=10 m','L=1 km','L=80 km','input')
grid on
subplot(1,2,2)
xlim([-5 5])
title('|E_t| con τ_0=1 ns numerica','FontSize',22)
xlabel('t [ns]','FontSize',16)      % plot dell'impulso d'ingresso numerico con tau_0= 1 ns
ylabel('|E_t|','FontSize',16)
grid on

figure(5)
subplot(1,2,1)
xlim([1.1599665e15 1.159978e15])
title('|E_ω| con τ_0=1 ns analitica','FontSize',22)
xlabel('ω [rad/s]','FontSize',16)        % plot dello spettro d'ingresso analitico con tau_0= 1 ns
ylabel('|E_ω|','FontSize',16)
legend('L=1 m','L=10 m','L=1 km','L=80 km','input')
grid on
subplot(1,2,2)
xlim([1.1599665e15 1.159978e15])
title('|E_ω| con τ_0=1 ns numerica','FontSize',22)
xlabel('ω [rad/s]','FontSize',16)        % plot dello spettro d'ingresso numerico con tau_0= 1 ns
ylabel('|E_ω|','FontSize',16)
grid on

figure(6)
subplot(1,2,1)
xlim([-0.7 0.7])
title('|E_t| con τ_0=10 ps analitica','FontSize',22)
xlabel('t [ns]','FontSize',16)      % plot dell'impulso d'uscita analitico con tau_0= 10 ps
ylabel('|E_t|','FontSize',16)
legend('L=1 m','L=10 m','L=1 km','L=80 km','input')
grid on
subplot(1,2,2)
xlim([-0.7 0.7])
title('|E_t| con τ_0=10 ps numerica','FontSize',22)
xlabel('t [ns]','FontSize',16)      % plot dell'impulso d'uscita numerico con tau_0= 10 ps
ylabel('|E_t|','FontSize',16)
grid on

figure(7)
subplot(1,2,1)
xlim([1.1596e15 1.16035e15])
title('|E_ω| con τ_0=10 ps analitica','FontSize',22)
xlabel('ω [rad/s]','FontSize',16)       % plot dello spettro d'uscita analitico con tau_0= 10 ps
ylabel('|E_ω|','FontSize',16)
legend('L=1 m','L=10 m','L=1 km','L=80 km','input')
grid on
subplot(1,2,2)
xlim([1.1596e15 1.16035e15])
title('|E_ω| con τ_0=10 ps numerica','FontSize',22)
xlabel('ω [rad/s]','FontSize',16)       % plot dello spettro d'uscita numerico con tau_0= 10 ps
ylabel('|E_ω|','FontSize',16)
grid on


%% parte 4

tau_0=[1e-9 5e-10 1e-10 1e-11]; % costanti di tempo [s]
alpha_dB=0.00023; % attenuazione della fibra [dB/m]
N=1e3; % fattore di campionamento
M=70; % fattore per l'ampiezza del segnale
L=linspace(0,8e4,N); % vettore lunghezze [m]
FWHM_tau=zeros(size(L)); % inizializzazione vettore di FWHM/tau_0

for k = 1:length(tau_0) % ciclo per gestire ogni costante di tempo
    t_end=M*tau_0(k); % limite per il vettore tempi
    dt=tau_0(k)/N; % tempo di campionamento
    t=(-t_end:dt:t_end); % vettore tempi [s]
    w=2*pi.*(linspace(-0.5/dt,0.5/dt,lung))+w_0; % vettore pulsazioni [rad/s]
    p_t=exp(-(t.^2/(2*tau_0(k)^2))); % impulso d'ingresso
    p_w=fftshift(fft(p_t))/lung; % spettro d'uscita
    for j = 1:length(L) % ciclo per gestire ogni lunghezza
        E_w=p_w.*exp(-0.5i.*Beta_0_secondo.*(w-w_0).^2.*L(j)).*10^(-alpha_db*L(j)/20); % spettro d'uscita
        E_t=ifft(E_w)*lung; % impulso d'uscita         
        half_amp=max(abs(E_t))/2; % mezza ampiezza del segnale 
        E_t_0=abs(E_t(1:(length(E_t)-1)/2)); % metà sinistra della gaussiana
        E_t_1=abs(E_t((length(E_t)-1)/2+1:end)); % metà destra della gaussiana
        [~, index_0] = min(abs(abs(E_t_0)-half_amp));  % indice di salita a mezza ampiezza
        [~, index_1_par] = min(abs(abs(E_t_1)-half_amp)); % indice di discesa a mezza ampiezza
        index_1=(length(E_t)-1)/2+index_1_par; % normalizzazione dell'indice di discesa
        FWHM_tau(j)=(t(index_1)-t(index_0))/tau_0(k); % valore del rapporto FWHM/tau_0
    end
    figure(8)
    hold on % plot FWHM
    plot(L,FWHM_tau,'LineWidth',0.5)
end

figure(8)
title('FWHM lungo la lunghezza','FontSize',22)
xlabel('L [m]','FontSize',16)              % plot FWHM 
ylabel('FWHM/τ_0','FontSize',16)
legend('τ_0=1 ns','τ_0=0.5 ns','τ_0=0.1 ns','τ_0=10 ps','Location','northwest')
grid on



%% parte 5

tau_0=1e-11; % costante di tempo [s]
L=[1000 80000]; % lunghezze [m]
N=1e7; % campioni
t=linspace(0,1e-3,N); % vettore tempi [s]

for j = 1:length(L) % ciclo per gestire ogni lunghezza
    chirp=-(Beta_0_secondo*L(j).*(t-Beta_0_primo*L(j)))./(tau_0^4+(Beta_0_secondo*L(j))^2); % chirp di frequenza [rad/s]
    figure(9)
    hold on % plot chirp
    plot(t,chirp,'LineWidth',0.5)
end

figure(9)
title('Chirp di frequenza con τ_0=10 ps','FontSize',22)
xlabel('t [s]','FontSize',16)         % plot del chirp di frequenza
ylabel('φ [rad/s]','FontSize',16)
legend('L=1 km','L=80 km','Location','east')
grid on




    