% roba=findobj(gca,'Type','Line');
% t=roba.XData;
% y_roba=roba.YData;

figure(1)
subplot(3,1,1)
plot(t,h0,'r-',t,h1,'b-',t,h2,'g-',t,h3,'m-','LineWidth',1.5)
title("P_{out} vs Time",FontSize=18)
xlabel("Time [ns]",FontSize=14)
ylabel("Output Power [mw]",FontSize=14)
xlim([9 20])
subplot(3,1,2)
plot(t,i0,'r-',t,i1,'b-',t,i2,'g-',t,i3,'m-','LineWidth',1.5)
title("Carrier Density vs Time",FontSize=18)
xlabel("Time [ns]",FontSize=14)
ylabel("Carrier Density [cm^{-3}]",FontSize=14)
xlim([9 20])
subplot(3,1,3)
plot(t,j0,'r-',t,j1,'b-',t,j2,'g-',t,j3,'m-','LineWidth',1.5)
title("\Delta\nu vs Time",FontSize=18)
xlabel("Time [ns]",FontSize=14)
ylabel("\Delta\nu [GHz]",FontSize=14)
xlim([9 20])
legend("I/I_{th}=2","I/I_{th}=1.05","I/I_{th}=6","I/I_{th}=10",FontSize=14)