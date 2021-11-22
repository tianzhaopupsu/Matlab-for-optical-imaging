close all

 destinwave=aa_interface_airfixbr5;
  destinwave1=aa_interface_waterfixbr5;
   destinwave2=aa_singleinsolumixbr;
 aa_airstats=zeros(10,10);
% 
% 
aa_airstats(1,1)=mean(destinwave(:,1));
aa_airstats(1,2)=std(destinwave(:,1));
aa_airstats(2,1)=mean(destinwave(:,2));
aa_airstats(2,2)=std(destinwave(:,2));
aa_airstats(3,1)=mean(destinwave(:,3));
aa_airstats(3,2)=std(destinwave(:,3));
aa_airstats(4,1)=mean(destinwave(:,4));
aa_airstats(4,2)=std(destinwave(:,4));
aa_airstats(5,1)=mean(destinwave(:,13));
aa_airstats(5,2)=std(destinwave(:,13));


aa_airstats(6,1)=mean(destinwave(:,5));
aa_airstats(6,2)=std(destinwave(:,5));
aa_airstats(7,1)=mean(destinwave(:,6));
aa_airstats(7,2)=std(destinwave(:,6));
aa_airstats(8,1)=mean(destinwave(:,7));
aa_airstats(8,2)=std(destinwave(:,7));
aa_airstats(9,1)=mean(destinwave(:,8));
aa_airstats(9,2)=std(destinwave(:,8));
aa_airstats(10,1)=mean(destinwave(:,9));
aa_airstats(10,2)=std(destinwave(:,9));
y = aa_airstats(1:5,1);
x = aa_airstats(6:10,1);
xneg = aa_airstats(6:10,2)/2;
xpos = aa_airstats(6:10,2)/2;
yneg = aa_airstats(1:5,2)/2;
ypos = aa_airstats(1:5,2)/2;


aa_airstats(1,3)=mean(destinwave1(:,1));
aa_airstats(1,4)=std(destinwave1(:,1));
aa_airstats(2,3)=mean(destinwave1(:,2));
aa_airstats(2,4)=std(destinwave1(:,2));
aa_airstats(3,3)=mean(destinwave1(:,3));
aa_airstats(3,4)=std(destinwave1(:,3));
aa_airstats(4,3)=mean(destinwave1(:,4));
aa_airstats(4,4)=std(destinwave1(:,4));
aa_airstats(5,3)=mean(destinwave1(:,13));
aa_airstats(5,4)=std(destinwave1(:,13));


aa_airstats(6,3)=mean(destinwave1(:,5));
aa_airstats(6,4)=std(destinwave1(:,5));
aa_airstats(7,3)=mean(destinwave1(:,6));
aa_airstats(7,4)=std(destinwave1(:,6));
aa_airstats(8,3)=mean(destinwave1(:,7));
aa_airstats(8,4)=std(destinwave1(:,7));
aa_airstats(9,3)=mean(destinwave1(:,8));
aa_airstats(9,4)=std(destinwave1(:,8));
aa_airstats(10,3)=mean(destinwave1(:,9));
aa_airstats(10,4)=std(destinwave1(:,9));
y1 = aa_airstats(1:5,3);
x1 = aa_airstats(6:10,3);
xneg1 = aa_airstats(6:10,4)/2;
xpos1 = aa_airstats(6:10,4)/2;
yneg1 = aa_airstats(1:5,4)/2;
ypos1 = aa_airstats(1:5,4)/2;

aa_airstats(1,5)=mean(destinwave2(:,1));
aa_airstats(1,6)=std(destinwave2(:,1));
aa_airstats(2,5)=mean(destinwave2(:,2));
aa_airstats(2,6)=std(destinwave2(:,2));
aa_airstats(3,5)=mean(destinwave2(:,3));
aa_airstats(3,6)=std(destinwave2(:,3));
aa_airstats(4,5)=mean(destinwave2(:,4));
aa_airstats(4,6)=std(destinwave2(:,4));
aa_airstats(5,5)=mean(destinwave2(:,13));
aa_airstats(5,6)=std(destinwave2(:,13));


aa_airstats(6,5)=mean(destinwave2(:,5));
aa_airstats(6,6)=std(destinwave2(:,5));
aa_airstats(7,5)=mean(destinwave2(:,6));
aa_airstats(7,6)=std(destinwave2(:,6));
aa_airstats(8,5)=mean(destinwave2(:,7));
aa_airstats(8,6)=std(destinwave2(:,7));
aa_airstats(9,5)=mean(destinwave2(:,8));
aa_airstats(9,6)=std(destinwave2(:,8));
aa_airstats(10,5)=mean(destinwave2(:,9));
aa_airstats(10,6)=std(destinwave2(:,9));
y2 = aa_airstats(1:5,5);
x2 = aa_airstats(6:10,5);
xneg2 = aa_airstats(6:10,6)/2;
xpos2 = aa_airstats(6:10,6)/2;
yneg2 = aa_airstats(1:5,6)/2;
ypos2 = aa_airstats(1:5,6)/2;



errorbar(x,y,yneg,ypos,xneg,xpos,'o','Color','red','MarkerFaceColor','red')
hold on

errorbar(x1,y1,yneg1,ypos1,xneg1,xpos1,'o','Color','green','MarkerFaceColor','green')
errorbar(x2,y2,yneg2,ypos2,xneg2,xpos2,'o','Color','blue','MarkerFaceColor','blue')
axis([-5 23 -0.1 0.8])
set(gca, 'Fontsize', 12)
ylabel('Weight of the component','Fontsize',15);
legend('Immoblized gQD in air','Immoblized gQD in solution','Diffusing gQD in solution','fontsize',12,'EdgeColor',[1 1 1]);
xlabel('Lifetime of the component (ns)','Fontsize',15);


