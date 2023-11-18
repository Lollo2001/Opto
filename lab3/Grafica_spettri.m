clc
clear all
close all

% grafica spettri
% imposta le correnti a cui sono stati misurati gli spettri e che si
% vogliono riportare sul grafico
Istart=1; % mA
Istop=40; % mA
Istep=1; %mA
Icurrent=Istart:Istep:Istop;
Mmeas=load('OSA_Temp20.000000.txt'); % carica il file contenente le misure di spettri ottici

num_punti=800; % numero di lunghezze d'onda a cui è stato misurato lo spettro ottico ad una corrente fissata
for icurr=1:length(Icurrent)
    figure(1)
    hold on
    plot(Mmeas((icurr-1)*num_punti+1:icurr*num_punti,1),Mmeas((icurr-1)*num_punti+1:icurr*num_punti,2))
    %ibias=icurr*Istep
    %pause
end
    

Istart=1; % mA
Istop=40; % mA
Istep=1; %mA
Icurrent=Istart:Istep:Istop;
Mmeas=load('OSA_Temp25.000000.txt'); % carica il file contenente le misure di spettri ottici

num_punti=800; % numero di lunghezze d'onda a cui è stato misurato lo spettro ottico ad una corrente fissata
for icurr=1:length(Icurrent)
    figure(2)
    hold on
    plot(Mmeas((icurr-1)*num_punti+1:icurr*num_punti,1),Mmeas((icurr-1)*num_punti+1:icurr*num_punti,2))
    %ibias=icurr*Istep
    %pause
end
    

Istart=1; % mA
Istop=40; % mA
Istep=1; %mA
Icurrent=Istart:Istep:Istop;
Mmeas=load('OSA_Temp30.000000.txt'); % carica il file contenente le misure di spettri ottici

num_punti=800; % numero di lunghezze d'onda a cui è stato misurato lo spettro ottico ad una corrente fissata
for icurr=1:length(Icurrent)
    figure(3)
    hold on
    plot(Mmeas((icurr-1)*num_punti+1:icurr*num_punti,1),Mmeas((icurr-1)*num_punti+1:icurr*num_punti,2))
    %ibias=icurr*Istep
    %pause
end
    

figure(1)
title("Spettri ottici con vari valori di corrente a T=20°C")
xlabel("λ [nm]")
ylabel("P [dBm]")

figure(2)
title("Spettri ottici con vari valori di corrente a T=25°C")
xlabel("λ [nm]")
ylabel("P [dBm]")

figure(3)
title("Spettri ottici con vari valori di corrente a T=30°C")
xlabel("λ [nm]")
ylabel("P [dBm]")