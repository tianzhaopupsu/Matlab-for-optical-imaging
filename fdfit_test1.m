
%cw is the center postion of the instrument respons function
%y is the fluoresence decay of the sample
%how many exponentials you want in this fitting, 1==1exp 2==2exp 3==3exp.
%setpsize is the time interval between pixels.
close all

a=61;
b=3786;
c=-2;

lb = [0 0 0 0 0 10 50 -10];
ub = [ 0.9 0.5 0.5 2 2 500 1000 100];
a0 = [0.5, 0.2, 0.01];            
tau0 = [2, 3, 50 450];
bkg = 0.01;

yo=aa_singleinsolu0916fix((a:b),2);


background=0;
%yo=yc(5:3897,:);
%yo=yc(53:3897,:);
timeaxis=aa_timeaxis(a:b,1);
funda=aa_irf_sol0907((a-c:b-c),1);
irf=funda;
stepsize=0.30517578125;
%cw=63;
cw=177;
%cw=23;
FWHM=0.9;%0.944  0.34 in ns
width=FWHM/stepsize;
fs=15; %font size
c=length(yo);
t=[0:(c-1)]';
y=yo./max(yo);
fitresult=zeros(5,5);

irf=irf./max(irf);
   


% fit the result using single expoential decay
t=timeaxis;

expfun = @(x,t,irf) model_trunc(t,irf,x(8),x(1:3),x(4:7)); 
fitresls =  @(x) ((((expfun(x,t,irf)/max(expfun(x,t,irf)))-y))./sqrt(y));



x0 = [a0, tau0, bkg];
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, x0,lb,ub);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, x0,x0,x0);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,[],[]);
ci = nlparci(x,resid,'jacobian',J);
chi2w=((expfun(x,t,irf)/max(expfun(x,t,irf)))-y).^2./y;
chi2a=max(yo).*sum(chi2w)./c;
t1=x(1)*x(4)+x(2)*x(5)+x(3)*x(6)+(1-x(1)-x(2)-x(3))*x(7)
% resultt=['t1 = ',num2str(x(2)) '+' num2str(x(2)-ci(2,1)) ' t2 = ',num2str(x(4)) '+' num2str(x(4)-ci(4,1))...
%     ' t3 = ',num2str(x(6)) '+' num2str(x(6)-ci(6,1))];
% sumpor = x(1)+x(3)+x(5)+x(7);
% resulta = ['t1 = ',num2str(x(1)/sumpor) '+' num2str((x(1)-ci(1,1))/sumpor) ' t2 = ',num2str(x(3)/sumpor) '+' num2str((x(3)-ci(3,1))/sumpor)...
%     ' t3 = ',num2str(x(5)/sumpor) '+' num2str((x(5)-ci(5,1))/sumpor)];
% fitresult(1,1)=x(2);
% fitresult(1,2)=x(2)-ci(2,1);
% fitresult(2,1)=x(4);
% fitresult(2,2)=x(4)-ci(4,1);
% fitresult(3,1)=x(6);
% fitresult(3,2)=x(6)-ci(6,1);
% fitresult(4,1)=x(8);
% fitresult(4,2)=x(8)-ci(8,1);
% fitresult(1,3)=x(1)/sumpor;
% fitresult(1,4)=(x(1)-ci(1,1))/sumpor;
% fitresult(2,3)=(x(3))/sumpor;
% fitresult(2,4)=(x(3)-ci(3,1))/sumpor;
% fitresult(3,3)=(x(5))/sumpor;
% fitresult(3,4)=(x(5)-ci(5,1))/sumpor;
% fitresult(4,3)=(x(7))/sumpor;
% fitresult(4,4)=(x(7)-ci(7,1))/sumpor;
% fitresult(5,1)=x(2)*x(1)/sumpor+x(4)*x(3)/sumpor+x(6)*x(5)/sumpor+x(8)*x(7)/sumpor;
% 

% plot (irs);
% hold on;
% txt = ['Chisq = ',num2str(B)];
% 
% text(140,1.2,txt,'VerticalAlignment','Cap');
% 
% plot(y);
% 
% plot(expfun(B,t)/max(expfun(B,t)))
% set(gca, 'YScale', 'log')
% axis([037 900 0.01 1.5])
% figure;
% plot(y);
% hold on
% 
% plot(expfun(B,t)/max(expfun(B,t)));
% set(gca, 'YScale', 'log');

figure;
pos1 = [0.1 0.43 0.85 0.5];
subplot('Position',pos1)
plot(timeaxis,y);
hold on
plot(timeaxis,expfun(x,t,irf)./max(expfun(x,t,irf)),'LineStyle','-','LineWidth',1.5);
plot(timeaxis,irf);
set(gca, 'Fontsize', fs-1)
ylabel('Normalized Distribution','Fontsize',fs);

xlabel('Time (ns)','Fontsize',fs);
 axis([0 1200 0.000006 1.2])
% text(140,1,resultt);
% text(140,0.6,resulta);
legend('Data','Fit','IRF','fontsize',14,'EdgeColor',[1 1 1]);
set(gca, 'YScale', 'log');
pos2 = [0.1 0.1 0.85 0.2];
subplot('Position',pos2)
plot(timeaxis,y-expfun(x,t,irf)./max(expfun(x,t,irf)),'LineStyle','-','LineWidth',1,'Marker','none',...
   'MarkerSize',0.5)
 axis([0 1200  -0.2 0.2])
 ylabel('Residuals','Fontsize',fs);

xlabel('Time (ns)','Fontsize',fs);
set(gca,'fontsize',fs-1)
hold off;

figure;
plot(timeaxis,y);
hold on
plot(timeaxis,expfun(x,t,irf)./max(expfun(x,t,irf)),'LineStyle','-','LineWidth',1.5);
plot(timeaxis,irf);
 axis([40 200 0.000006 1.2])
 set(gca, 'YScale', 'log');
 hold off;
expfunr=expfun(x,t,irf);
aexpfunrn=expfun(x,t,irf)/max(expfun(x,t,irf));


function y = model_trunc( t, irf, bkg,  a, tau )
% MODEL_TRUNC generates lifetime model, truncated to valid-chanel range.
%  t, delay time array, parameter
%  irf, instrument response function array, parameter
%  bkg, background estimated from t<0, parameter
%  t0, t-zero shift, variable
%  a, amplitude array, variable
%  tau, lifetime array, variable
  tt = t;
 
  y0 = 0;
  aa = [a, 1-sum(a)];
  for i=1:length(tau)
    y0 = y0 + aa(i).*exp(-tt./tau(i));
  end
y0=y0;
  y1 = conv( y0, irf )+bkg;
  y1 = y1 ./ max(y1);
  y = y1(1:length(t));
end

