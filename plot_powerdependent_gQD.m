close all
fs=15;
plot(tax,testair(:,1)/max(testair(:,1)))
hold on
plot(tax,testair(:,2)/max(testair(:,2)))
plot(tax,testair(:,3)/max(testair(:,3)))
plot(tax,testair(:,4)/max(testair(:,4)))
plot(tax,testair(:,5)/max(testair(:,5)))
plot(tax,testair(:,6)/max(testair(:,6)))
plot(tax,testair(:,7)/max(testair(:,7)))
plot(tax,testair(:,8)/max(testair(:,8)))
plot(tax,testair(:,9)/max(testair(:,9)))
plot(tax,testair(:,10)/max(testair(:,10)))
plot(tax,testair(:,11)/max(testair(:,11)),'color','black')
set(gca, 'YScale', 'log');
legend('Bead 1','Bead 2','Bead 3','Bead 4','Bead 5','Bead 6','Bead 7','Bead 8','Bead 9','Bead 10','Bead ensemble',...
    'fontsize',10,'EdgeColor',[1 1 1]);
yl = ylim;
set(gca, 'Fontsize', fs-1)
ylabel('Counts','Fontsize',fs);
xlabel('Time (ns)','Fontsize',fs);
axis([0 150 1 yl(2)])


