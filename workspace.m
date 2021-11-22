close all
% % % %b=0;
% % % % for i=1:202
% % % %     
% % % %     if batchresult_water(i,7)>0
% % % %         
% % % %         b=b+1;
% % % %     end
% % % %     
% % % %     
% % % %     
% % % % end
% % % mea=mean(aa_singleinsolusumintHaw)
% % % st=std(aa_singleinsolusumintHaw)
% % % destwave=aa_diameter0916(1:233,1);
% % % 
% % % fprintf( 1, ['avg' num2str(mean(destwave(:,1))) 'std' num2str(std(destwave(:,1)))] );
% % % 
% % % 
% % % 
% % % % timeaxis=0:0.1220703125:500;
% % % % timeaxis=timeaxis';
% % % %  i=5;
% % % %  plot(bb(:,2),bb(:,1));
% % % % % 
% % % % % for i=1:202
% % % %  timec(:,1)=timeaxisf;
% % % %  timec(:,2)=Surfacesolution;
% % 
% % % % %ap='C:\Users\Chopin Pro\Desktop\Lifetime\Interface_water_08172021a\';
% % % % % ap='IRF_for_bulk';
% % % % % 
% % % % % bp='.txt';
% % % % % fName=append(ap,bp);
% % % %  save('Bulk_surface_in_solution.txt','timec','-ascii');
% % % % % end
% % %  plot(aa_IRF125_SPC,'Color','black');
% fs=15;
%  plot(aa_timeaxis,aa_bulkinsol0907);
%  hold on
%   plot(aa_timeaxis,aa_800*15);
%   plot(aa_timeaxis-16,aa_irf(:,1));
%   plot(aa_timeaxis,aa_irf_sol0907(:,1));
% ylabel('Counts','Fontsize',fs);
%  xlabel('Time (ns)','Fontsize',fs);
% set(gca, 'YScale', 'log');
% % % % for i=2:165
% % % %     %plot(aa_timeaxis,aa_interface_airfix(:,i));
% % % %     plot(aa_timeaxis,aa_singleinsolu0916fix(:,i)./max(aa_singleinsolu0916fix(:,i)));
% % % % end
% % %plot(timeax125_SPC,Rd6g125_SPC,'Color','red','Linewidth',2);
% % % plot(aa_interface_waterfixbr5(:,11),aa_interface_waterfixbr5(:,12),'LineStyle','none','Marker','*','MarkerSize',3,'MarkerEdgeColor','green');
% % % % % % plot(IRFlaser);
% % % % % % plot(Interface_air(:,1));
% % % % % % plot(IRF_lasera);
% % % % % % % 
% % % plot(FWHMMaxis/23.55,100-(100*(resultss(:,1)-3.85)/3.85))
% 
% % % % % % % % 
% % % % %
% % % xlabel('FWHM (ns)','Fontsize',fs);
% axis([0 220 0 70000])
% % % % % % % % text(140,1,resultt);
% % % % % % % % text(140,0.6,resulta);
% legend('400 nm excitation','800 nm excitation','IRF 400nm','IRF 800nm','fontsize',14,'EdgeColor',[1 1 1]);
% % % set(gca, 'Fontsize', 14)
% % %  %
% % % % axis([0 4  0 10])
% % % axis([10 1110  10 100000])
% % % 
% % fs=15;
% %  files=dir('\\192.168.1.101\newdata\tian\Y2021\M11\D01\Zoff\');
% % % %for i=3:length(files)
% % for i=3:550
% % close all
% % a=['\\192.168.1.101\newdata\tian\Y2021\M11\D01\Zoff\' files(i).name];
% %  fa=importdata(a);
% % x=[0 35]; 
% % y=[0 35];
% % imagesc(x,y,fa)
% % caxis([0 20]);
% % set(gca,'fontsize',15);
% %  xlabel('X (\mum)','FontSize',15)
% % ylabel('Y (\mum)','FontSize',15)
% % c = colorbar;
% %     
% %        set(c,'fontsize',fs-1);
% % c.Label.String = ['Photon counts (kcps)'];
% %  saveas(gca,num2str(i),'png');
% %   b=i
% %   end  
% % % a='\\192.168.1.101\newdata\tian\Y2021\M11\D01\Mean_intensity21-11-01_005924';
% % % fid=fopen(a,'r','b');
% % % chan1=fread(fid,[2,inf],'double')';
% % 
% % 
% % 
% % % aaa=chan1(1:5000,1:3);
% % % plot(chan1(:,2))
% % % aa_IRF125_SPC=fa;
% % % % % % % aaa5000(479:980,3)=aaa5000(479:980,2)/3;
% % % % % % % %t=1:0.305:1250;
% % % plot(zvar(:,3),zvar(:,1)-zvar(1,1),'LineWidth',1);
% % % hold on
% % % plot(zvar(:,3),zvar(:,2)-zvar(2,2),'LineWidth',1);
% % % 
% % % % % plot(a4,'LineStyle','-','LineWidth',1);
% % % legend('Z correction on','Z correction off',...
% % %  'fontsize',14,'EdgeColor',[1 1 1]);
% % % % 
% % % % % % % % 
% % % ylabel('Z position (\mum)','Fontsize',15) 
% % % xlabel('Time (min)','Fontsize',15)

% % % % % sump(1,:)=sum(aaa5000(98:122,1));
% % % % % sump(2,:)=sum(aaa5000(489:980,3));
% % % % % sump(3,:)=sum(aaatest(610:3069,3));
% % % % 
% % % 
% % % % % 
% % % % % % % % % 
%fa=importdata('\\192.168.1.101\newdata\tian\Y2021\M11\D19\FLIM_191121_164543\Lifetime21-11-19_16  1 53');

plot(taxis,interlf);
% hold on
% %plot(aa_t,aa_QD_SPCM);
%set(gca, 'YScale', 'log')
set(gca,'fontsize',12);
 ylabel('Counts','Fontsize',15) 
 xlabel('Time (ns)','Fontsize',15)
%axis([0 450 0 8000])
% legend('LabView Measured','SPCM software','fontsize',14,'EdgeColor',[1 1 1]);
% % % % % % % % % % % BW=a>6;
% % % % % % % % % % % centers=imfindcircles(BW,[1 20]);
%  x = [0 40];
% y = [0 40];
% imagesc(x,y,FLIMimag);
% %imagesc(x,y,fig5);
%  caxis([10 50]);
% colorbar;
%  c = colorbar;
%  c.Label.String = 'Lifetime (ns)';
% set(gca,'fontsize',15);
%   xlabel('X (\mum)','FontSize',15)
%         ylabel('Y (\mum)','FontSize',15)

% % % % plot(aa_singleinsolu0907fixbr5(:,11))
% % % % hold on 
% % % %plot(aexpfunrn1)
