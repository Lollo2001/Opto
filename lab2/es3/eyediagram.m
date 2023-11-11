clc

%{
Genera il grafico, seleziona la linea, vai su command windows e manda questi comandi
h=gco
time=get(h, 'xdata')
Py=get(h, 'ydata')

Oppure per il grafio a 0.1Ghz importa i vettori time.mat e Py.mat
%}


t0=5; % ns
T_bit=10; % ns

figure()
for i=1:64
    itime=find(time>=t0-T_bit/2+2*T_bit*(i-1) & time<=t0+T_bit*3/2+2*T_bit*(i-1));  % vettore con indici dei punti che stanno dentro alla finestra scelta 
    time1=time(itime)-2*T_bit*(i-1);    % sposto la base tempi della finestra qualsiasi sulla base tempi della prima finestra
    Pf1=Py(itime);  % valori della potenza della finestra scelta
    
    hold on
    plot(time1,Pf1, 'b')
end
