
figure
subplot(3,1,1)
plot(x1,y1,'r','LineWidth',2)
hold on
plot(x3,y3,'LineWidth',2,'Color',[0.07 0.62 1.00])
hold on
plot(x5,y5,'LineWidth',2,'Color',[1.00 0.00 1.00])
hold on
plot(x7,y7,'LineWidth',2,'Color',[0.00 1.00 0.00])

subplot(3,1,2)
plot(x2,y2,'r','LineWidth',2)
hold on
plot(x4,y4,'LineWidth',2,'Color',[0.07 0.62 1.00])
hold on
plot(x6,y6,'LineWidth',2,'Color',[1.00 0.00 1.00])
hold on
plot(x8,y8,'LineWidth',2,'Color',[0.00 1.00 0.00])

subplot(3,1,1)
xlim([1 30])
title('P_{out} against Time','FontSize',11)
xlabel('Time [ns]','FontSize',11)    
ylabel('Output Power [mw]','FontSize',11)
grid on

subplot(3,1,2)
xlim([1 30])
title('Carrier Density against Time','FontSize',11)
xlabel('Time [ns]','FontSize',11)     
ylabel('Carrier Density [cm^{-}3]','FontSize',11)
grid on

legend(' Caso in esame',' Caso: 2*C',' Caso: 3*C',' Caso: 4*C', 'Location','south','FontSize',14)