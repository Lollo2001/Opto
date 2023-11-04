clc
clear all
close all

%% 1 - Valutazione n_eff(d)
sample = 50;

n_cladding=1.44;
n_core=3.48;
lambda=1.55*10^-6;
omega=2*pi*3e8/lambda;
mu=4*pi*10^-7;
m0=0;
m1=1;
m2=2;
k_0=2*pi/lambda;

n_eff=linspace(n_cladding,n_core-0.0001,sample);
gamma= k_0*(sqrt(n_eff.^2 -n_cladding^2));

k_x= k_0*(sqrt(n_core^2-n_eff.^2));

d0=(2*atan(gamma./k_x)+m0*pi)./k_x;
d1=(2*atan(gamma./k_x)+m1*pi)./k_x;
d2=(2*atan(gamma./k_x)+m2*pi)./k_x;

figure (1)
plot( d0, n_eff,'m--', d1, n_eff, 'g-', d2, n_eff, 'b')
xlabel('d, Thickness (m)');
ylabel('n_eff, effective refractive index ');
legend('m0','m1','m2');
title('n_eff(d)');
xlim([0, 3e-6]);
ylim([n_cladding, n_core]);
%% punto 2

d=4.5e-7; %value for having only te0 and te1 modes
[~, index_0] = min(abs(d0-d));  % trovo l'indice del dato del vettore d0 con il valore più vicino a d
[~, index_1] = min(abs(d1-d));  % trovo l'indice del dato del vettore d1 con il valore più vicino a d
n_eff_0 = n_eff(index_0);      % valore efficace modo 0
n_eff_1 = n_eff(index_1);      % valore efficace modo 1
beta_0 = k_0*n_eff_0;       % beta efficace modo 0
beta_1 = k_0*n_eff_1;       % beta efficace modo 1
k_x_0 = k_0*(sqrt(n_core^2-n_eff_0^2));
k_x_1 = k_0*(sqrt(n_core^2-n_eff_1^2));
gamma_0 = k_0*(sqrt(n_eff_0^2 -n_cladding^2));
gamma_1 = k_0*(sqrt(n_eff_1^2 -n_cladding^2));

% Ey
syms x;
Ey_0 = piecewise( x>0, exp(-gamma_0.*x), -d<=x<=0, (cos(k_x_0.*x) - (gamma_0/k_x_0)*sin(k_x_0.*x)), x<-d, (cos(k_x_0.*(-d)) - (gamma_0/k_x_0)*sin(k_x_0.*(-d)))*exp(gamma_0*(x+d)) );
figure(2)
hold on
grid on
fplot(Ey_0, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("E_y modo 0")

Ey_1 = piecewise( x>0, exp(-gamma_1.*x), -d<=x<=0, (cos(k_x_1.*x) - (gamma_1/k_x_1)*sin(k_x_1.*x)), x<-d, (cos(k_x_1.*(-d)) - (gamma_1/k_x_1)*sin(k_x_1.*(-d)))*exp(gamma_1*(x+d)) );
figure(3)
hold on
grid on
fplot(Ey_1, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("E_y modo 1")

% Hz    [ dEy/dx = -j*omega*mu*Hz ]
Hz_0 = diff(Ey_0, x)/(omega*mu);
figure(4)
hold on
grid on
fplot(Hz_0, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_z modo 0")

Hz_1 = diff(Ey_1, x)/(omega*mu);
figure(5)
hold on
grid on
fplot(Hz_1, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_z modo 1")

% Hx    [ -j*beta*Ey = j*omega*mu*Hx ]
Hx_0 = -(beta_0*Ey_0)/(omega*mu);
figure(6)
hold on
grid on
fplot(Hx_0, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_x modo 0")

Hx_1 = -(beta_1*Ey_1)/(omega*mu);
figure(7)
hold on
grid on
fplot(Hx_1, [-2*d d])
xline(-d, '--k');
xline(0, '--k');
title("H_x modo 1")
%{


%% punto 3
% ha senso utilizzare un range di d0 che va da 0 fino al massimo di 2/3e-6

Gamma = ones(1,length(n_eff));

for i = 1:length(n_eff)     % prendo uno a uno tutti i valori di d così so il corrispettivo n_eff e posso calcolarmi Ey
    syms x;
    Ey = piecewise( x>0, exp(-gamma(i).*x), -d0(i)<=x<=0, (cos(k_x(i).*x) - (gamma(i)/k_x(i))*sin(k_x(i).*x)), x<-d0(i), (cos(k_x(i).*(-d0(i))) - (gamma(i)/k_x(i))*sin(k_x(i).*(-d0(i))))*exp(gamma(i)*(x+d0(i))) );
    Gamma(i) = int(Ey^2, x, -d0(i), 0)/int(Ey^2, x, -inf, inf);    
end 

figure(8)
plot(d0, Gamma)
title("{Γ_TE0} in funzione dello spesso del core");
xlabel("thickness");
ylabel("Γ_TE0");



%% punto 4

alpha_core = 0.3;
alpha_cladding = 0.15;

alpha_modal = Gamma.*alpha_core + (1-Gamma).*alpha_cladding;    % alpha in cm^-1
alpha_modal_db = alpha_modal.*10*log10(exp(1));                 % alpha in dB/cm
figure(9)
plot(d0, alpha_modal)
figure(10)
semilogx(d0, alpha_modal_db)    % da capire se è meglio avere uno/entrambi gli assi logaritmici oppure no

%slide 13 pacco 2 reflection coefficient=d/we che minchia è we?

%% punto 5

len = 0.5; % guida lunga 0.5 cm
insertion_loss = alpha_modal_db*len;    % non sono sicuro ma da quello che ho capito la formula è -10*log(P_out/P_in)
figure(11)
plot(d0, insertion_loss)

%% punto 6
[~, index_Gamma] = min(abs(Gamma-0.75));  % trovo l'indice del dato del vettore Gamma con il valore più vicino a 0.75
thickness = d0(index_Gamma);

%% punto 7
lambda_7 = linspace(850, 1620, 1e4);    % vettore lambda da 850 nm a 1620 nm


%}















