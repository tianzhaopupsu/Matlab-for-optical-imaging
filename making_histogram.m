close all
destinationwave=aa_interface_waterfixint;
destinationwave1=aa_interface_airfixHaw;
figure;

[m,n]=size(destinationwave);
samplenumber=m;
fs=15;
 binsize=round(1+3.32*log10(samplenumber));
 histogram(destinationwave,binsize,'FaceColor','red')

axis([32 35 0 100])
ylabel('Occurrence','Fontsize',fs);
xlabel('Total lifetime (ns)','Fontsize',fs);
%xlabel('\chi^2','Fontsize',fs);
set(gca, 'Fontsize', fs-1)
%set(gca, 'YScale', 'log')

figure;

[m,n]=size(destinationwave1);
samplenumber=m;
fs=15;
 binsize=round(1+3.32*log10(samplenumber));
 histogram(destinationwave1,binsize,'FaceColor','red')

axis([1 200 0 100])
ylabel('Occurrence','Fontsize',fs);
%xlabel('Total lifetime (ns)','Fontsize',fs);
xlabel('Lifetime (ns)','Fontsize',fs);
set(gca, 'Fontsize', fs-1)
%set(gca, 'YScale', 'log')