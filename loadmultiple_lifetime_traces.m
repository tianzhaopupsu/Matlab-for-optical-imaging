close all


for i = 1:165

filenumber=string(i);


%ap='\\192.168.1.101\newdata\tian\Y2021\M09\D13\Interface_water_09102021\D0909_solids_sp';
ap='C:\Users\Chopin Pro\Desktop\Lifetime\D09142021_s35_p';

bp='.asc';
fName=append(ap,filenumber,bp);
fa=importdata(fName);
% fa=importdata('\\192.168.1.101\newdata\tian\Y2021\M08\D17\Water_rasterscan\NormalizedImages_180821_023206\Z Pixel 0\APD X Left.csv');
% imagesc(fa)
% x = [0 10];
% y = [0 10];
% imagesc(x,y,fa);
% caxis([0 30]);
% colorbar;
% c = colorbar;
% c.Label.String = 'Counts (kcps)';
% set(gca,'fontsize',15);
%   xlabel('X (\mum)','FontSize',15)
%         ylabel('Y (\mum)','FontSize',15)
%testwater(:,a)=fa(:,2);
aa_test(:,i)=fa.data;
end
% %b=sum(testwater-Interface_water);
% %a=sum(b);