clc
close all
clear

% question 1
k_q_0=0.01;
k_q_1=0.1;
n_eff=2.6;
n_g=3.2;
loss=50; %db/m
lambda=1.55*10^(-6);
c=3*10^(8); %m/s
r_r=8*10^(-6); %m

lambda_v=linspace(1.50*10^(-6),1.60*10^(-6),10^5);
f=c./lambda_v;
w=2*pi.*f;

%per le due gamme diverse, calcolo del coefficiente per il campo elettrico
%che esce dalla guida dritta, non entra nell'anello
t_0=sqrt(1-k_q_0); 
t_1=sqrt(1-k_q_1); 

L=r_r*2*pi; %lunghezza dell'anello
aL=5*10^(-6);
a=exp(aL/2); %perdite
beta=w.*n_eff/c; %coefficiente di propagazione

%per le due gamme diverse, calcolo del coefficiente di transmissione
T_thru_0=(t_0-t_0*a.*exp(beta.*1i*L))./(1-a*t_0^2.*exp(beta.*1i*L)); %coefficiente di trasmissione
T_drop_0=(-k_q_0*sqrt(a).*sqrt(exp(beta.*1i*L)))./(1-a*t_0^2.*exp(beta.*1i*L)); %coefficiente di trasmissione

T_thru_1=(t_1-t_1*a.*exp(beta.*1i*L))./(1-a*t_1^2.*exp(beta.*1i*L)); %coefficiente di trasmissione
T_drop_1=(-k_q_1*sqrt(a).*sqrt(exp(beta.*1i*L)))./(1-a*t_1^2.*exp(beta.*1i*L)); %coefficiente di trasmissione

%power transmission spectrum at the through and drop
T0T=abs(T_thru_0).^2;
T0D=abs(T_drop_0).^2;
T1T=abs(T_thru_1).^2;
T1D=abs(T_drop_1).^2;

figure (1)
subplot(2,1,1)
plot(lambda_v,T0T)
title("Power at the through")
xlabel("Lambda [m]")
ylabel("Power")
grid on
subplot(2,1,2)
plot(lambda_v,T0D)
title("Power at the drop")
xlabel("Lambda [m]")
ylabel("Power")
grid on

figure (2)
subplot(2,1,1)
plot(lambda_v,T1T)
title("Power at the through")
xlabel("Lambda [m]")
ylabel("Power")
grid on
subplot(2,1,2)
plot(lambda_v,T1D)
title("Power at the drop")
xlabel("Lambda [m]")
ylabel("Power")
grid on

% we can use the ring as a sensor when k^2=0.1 because it's more sensible

%%
% question 2

lambda_m=lambda;
alfa_analyte=linspace(0,3000,10^6);
Gamma_cl_0=0.01;
Gamma_cl_1=0.05;
Gamma_cl_2=0.1;

w_FSR=2*pi*c/(n_g*L);

%perdita dal dB/m al 1/m
loss_m_q=10^((loss*L)/10);
loss_m=log(loss_m_q)/L; 

%fattore di qualità
Q_0=(pi*n_g*L*sqrt(t_0^2.*exp(-(loss_m+alfa_analyte*Gamma_cl_0).*L/2)))./((1-t_0^2.*exp(-(loss_m+alfa_analyte*Gamma_cl_0)*L/2))*lambda_m);
Q_1=(pi*n_g*L*sqrt(t_0^2.*exp(-(loss_m+alfa_analyte*Gamma_cl_1).*L/2)))./((1-t_0^2.*exp(-(loss_m+alfa_analyte*Gamma_cl_1)*L/2))*lambda_m);
Q_2=(pi*n_g*L*sqrt(t_0^2.*exp(-(loss_m+alfa_analyte*Gamma_cl_2).*L/2)))./((1-t_0^2.*exp(-(loss_m+alfa_analyte*Gamma_cl_2)*L/2))*lambda_m);

figure (3)
plot(alfa_analyte,Q_0,alfa_analyte,Q_1,alfa_analyte,Q_2)
title("Quality factor")
legend("Q0","Q1","Q2","Location","south")
xlabel("alfa analyte [m^{-1}]")
ylabel("Q")
grid on

%%
% question 3
%calcolo sensibilità per ogni gamma e limite di rilevamento
Gamma_cl_v=linspace(0.01,0.1,1e4);
S=Gamma_cl_v.*lambda./n_eff;

S_0=Gamma_cl_0*lambda/n_eff;
DL_0=lambda./(Q_0.*S_0);

S_1=Gamma_cl_1*lambda/n_eff;
DL_1=lambda./(Q_1.*S_1);

S_2=Gamma_cl_2*lambda/n_eff;
DL_2=lambda./(Q_2.*S_2);

figure (4)
plot(Gamma_cl_v,S,"k")
xlabel("\Gamma_c_l_a_d_d_i_n_g")
ylabel("S [m]")
grid on
hold on
plot(Gamma_cl_0,S_0,"*",Gamma_cl_1,S_1,"*",Gamma_cl_2,S_2,"*")
legend("","\Gamma_c_l_a_d_d_i_n_g=0.01","\Gamma_c_l_a_d_d_i_n_g=0.05","\Gamma_c_l_a_d_d_i_n_g=0.1","Location","southeast")
title("Sensitivity")
grid on

figure (5)
plot(alfa_analyte,DL_0,alfa_analyte,DL_1,alfa_analyte,DL_2)
title("Detection Limit")
legend("DL0","DL1","DL2","location","east")
xlabel("alfa analyte [m^{-1}]")
ylabel("DL")
grid on


