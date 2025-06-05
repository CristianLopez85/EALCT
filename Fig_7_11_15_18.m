%% Fig. 7
figure; set(gcf,'Position',[20 100 600 300]);

plot(Res(:,1),Res(:,2),'k','Marker','diamond','MarkerFaceColor','k')
hold on;
plot(Res(:,1),Res(:,3),'Color',[0 0.4 0],'Marker','s','MarkerFaceColor',[0 0.4 0])
plot(Res(:,1),Res(:,4),'b','Marker','o','MarkerFaceColor','b')
plot(Res(:,1),Res(:,5),'m','Marker','>','MarkerFaceColor','m')
plot(Res(:,1),Res(:,6),'r','Marker','^','MarkerFaceColor','r')

xlim([-4.2 15.2]); 
set(gca, 'YScale', 'log')
xlabel('SNR (dB)'); ylabel('Error of signals estimation');
ylim([0.03 1.3])
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
legend('STFT + R-I','GLCT + R-I','ALCT + R-I','ALCT','EALCT','Location','best','Fontsize',10)
%% Fig. 11
figure; set(gcf,'Position',[20 100 600 300]);

plot(Res(:,1),Res(:,2),'k','Marker','diamond','MarkerFaceColor','k')
hold on;
plot(Res(:,1),Res(:,3),'Color',[0 0.4 0],'Marker','s','MarkerFaceColor',[0 0.4 0])
plot(Res(:,1),Res(:,4),'b','Marker','o','MarkerFaceColor','b')
plot(Res(:,1),Res(:,5),'m','Marker','>','MarkerFaceColor','m')
plot(Res(:,1),Res(:,6),'r','Marker','^','MarkerFaceColor','r')

xlim([-4.2 15.2]); 
set(gca, 'YScale', 'log')
xlabel('SNR (dB)'); ylabel('Error of signals estimation');
ylim([0.025 1.2])
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
legend('STFT + R-I','GLCT + R-I','ALCT + R-I','ALCT','EALCT','Location','best','Fontsize',10)
%% Fig. 15
figure; set(gcf,'Position',[20 100 600 300]);

plot(Res(:,1),Res(:,2),'k','Marker','diamond','MarkerFaceColor','k')
hold on;
plot(Res(:,1),Res(:,3),'Color',[0 0.4 0],'Marker','s','MarkerFaceColor',[0 0.4 0])
plot(Res(:,1),Res(:,4),'b','Marker','o','MarkerFaceColor','b')
plot(Res(:,1),Res(:,5),'m','Marker','>','MarkerFaceColor','m')
plot(Res(:,1),Res(:,6),'r','Marker','^','MarkerFaceColor','r')

xlim([-4.2 15.2]); 
set(gca, 'YScale', 'log')
xlabel('SNR (dB)'); ylabel('Error of signals estimation');
ylim([0.03 2])
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
legend('STFT + R-I','GLCT + R-I','ALCT + R-I','ALCT','EALCT','Location','best','Fontsize',10)

%% Fig. 18
figure; set(gcf,'Position',[20 100 600 300]);

plot(Res(:,1),Res(:,2),'k','Marker','diamond','MarkerFaceColor','k')
hold on;
plot(Res(:,1),Res(:,3),'Color',[0 0.4 0],'Marker','s','MarkerFaceColor',[0 0.4 0])
plot(Res(:,1),Res(:,4),'b','Marker','o','MarkerFaceColor','b')
plot(Res(:,1),Res(:,5),'m','Marker','>','MarkerFaceColor','m')
plot(Res(:,1),Res(:,6),'r','Marker','^','MarkerFaceColor','r')

xlim([-4.2 15.2]); 
set(gca, 'YScale', 'log')
xlabel('SNR (dB)'); ylabel('Error of signals estimation');
ylim([0.05 1.1])
set(findall(gcf,'-property','FontSize'),'FontSize',12, 'FontName', 'Times New Roman')
legend('STFT + R-I','GLCT + R-I','ALCT + R-I','ALCT','EALCT','Location','best','Fontsize',10)
