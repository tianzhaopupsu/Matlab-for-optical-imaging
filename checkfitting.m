close all
number=37;
max=70000;
irfsca=20;
a=61;
b=3786;
checklifetime=aa_singleinsolu0907fix(a:b,number);
checkfit=aa_singleinsolu0907fixbfr5(:,number);
checkirf=aa_irf_sol0907(a:b,1);
timeaxis=aa_timeaxis(a:b,1);
fs=15;
h=figure;

plot(timeaxis,checklifetime);
hold on
plot(timeaxis,checkfit,'LineStyle','-','LineWidth',1.5);
plot(timeaxis,checkirf./irfsca);
set(gca, 'Fontsize', fs-1)
ylabel('Normalized Distribution','Fontsize',fs);

 axis([0 1000 1 max])
% text(140,1,resultt);
% text(140,0.6,resulta);
legend('Data','Fit','IRF','fontsize',14,'EdgeColor',[1 1 1]);
set(gca, 'YScale', 'log');

hold off;
saveas(h,['e' num2str(number) '.jpg'])

g=figure;


plot(timeaxis,checklifetime);
hold on
plot(timeaxis,checkfit,'LineStyle','-','LineWidth',1.5);
plot(timeaxis,checkirf./irfsca);
set(gca, 'Fontsize', fs-1)
ylabel('Normalized Distribution','Fontsize',fs);

 axis([35 200 1 max])
% text(140,1,resultt);
% text(140,0.6,resulta);
legend('Data','Fit','IRF','fontsize',14,'EdgeColor',[1 1 1]);
set(gca, 'YScale', 'log');
set(gca,'fontsize',fs-1)
hold off;
saveas(g,['e' num2str(number)  'a' '.jpg'])
chi2=aa_singleinsolu0907fixbr5(number,11)
