%% caso T = 20

clear 
close all 
clc


filep=load('Parte2_dir1\PI_Temp20,00_IStart0,10_IStop12,00_IStep0,10.txt');
filev=load('Parte2_dir1\VI_Temp20,00_IStart0,10_IStop12,00_IStep0,10.txt');
figure(1)
subplot(1,2,1)
plot(filep(:,1),10.*filep(:,2))
subplot(1,2,2)
plot(filev(:,1),filev(:,2))

figure(1)
subplot(1,2,1)
title("Potenza d'uscita")
xlabel("I [mA]")
ylabel("P [W]")
subplot(1,2,2)
title("Tensione d'uscita")
xlabel("I [mA]")
ylabel("V [V]")

%% caso no T

filep=load('Parte2_dir2\PI_IStart0,00_IStop50,00_IStep1,00.txt');
filev=load('Parte2_dir2\VI_IStart0,00_IStop50,00_IStep1,00.txt');
figure(2)
subplot(1,2,1)
plot(filep(:,1),10.*filep(:,2))
subplot(1,2,2)
plot(filev(:,1),filev(:,2))

figure(2)
subplot(1,2,1)
title("Potenza d'uscita")
xlabel("I [mA]")
ylabel("P [W]")
subplot(1,2,2)
title("Tensione d'uscita")
xlabel("I [mA]")
ylabel("V [V]")

%% caso vari T 

for k = 20:30
    filep=load("Parte2_dir3\PI_Temp"+k+",00_IStart0,00_IStop40,00_IStep0,50.txt");
    filev=load("Parte2_dir3\VI_Temp"+k+",00_IStart0,00_IStop40,00_IStep0,50.txt");
    figure(3)
    hold on
    subplot(1,2,1)
    plot(filep(:,1),10.*filep(:,2))
    figure(3)
    hold on
    subplot(1,2,2)
    plot(filev(:,1),filev(:,2))
end

figure(3)
subplot(1,2,1)
title("Potenza d'uscita")
xlabel("I [mA]")
ylabel("P [W]")
subplot(1,2,2)
title("Tensione d'uscita")
xlabel("I [mA]")
ylabel("V [V]")
