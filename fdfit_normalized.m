
%cw is the center postion of the instrument respons function
%y is the fluoresence decay of the sample
%how many exponentials you want in this fitting, 1==1exp 2==2exp 3==3exp.
%setpsize is the time interval between pixels.
close all

a=57;
b=250;
cc=2;
yo=binedlf(a:b);
timeaxis=taxisbin(a:b,1);
funda=taxis(a-cc:b-cc,1);
expn=3;
irff=2;
background=0;
%yo=yc(5:3897,:);
%yo=yc(53:3897,:);

irf=funda;
stepsize=timeaxis(2)-timeaxis(1);
%cw=63;
cw=25;
%cw=23;
FWHM=0.8;%0.944  0.34 in ns
width=FWHM/stepsize;
fs=15; %font size
c=length(yo);
t=[0:(c-1)]';
y=yo./max(yo);
y(y==0)=1;
fitresult=zeros(4,4);

if irff==1
irf=irf./sum(irf);
    irs=irf;

end

if irff==2
irs=exp(-((t-cw)/width).^2);
end

% fit the result using single expoential decay
if expn==1

expfun = @(tau,t) (real(ifft((fft(irs)*ones(1,size(y,2))).*fft(exp(-(t/(tau(1)/stepsize)))+tau(2)))));


fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));

start_point =[3.8 background]; 
start_point1 =[3.6 background*0]; 
start_point2 =[3.9 background*150000]; 

%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point, start_point1,start_point2);
ci = nlparci(x,resid,'jacobian',J);
chi2w=((expfun(x,t)/max(expfun(x,t)))-y).^2./y;
chi2=max(yo).*sum(chi2w)./c;
resultt=['t = ',num2str(x(1)) '+' num2str(x(1)-ci(1,1))];
resulta = num2str(1);
fitresult(1,1)=x(1);
fitresult(1,2)=x(1)-ci(1,1);
end 

% fit the result using double expoential decay
if expn==2

expfun = @(a,t) (real(ifft((fft(irs)*ones(1,size(y,2))).*fft(a(1)*exp(-(t/(a(2)/stepsize)))+(1-a(1))*exp(-(t/(a(3)/stepsize)))+a(4)))));

fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));
start_point =[1,35,600,background]; 
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,[],[]);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point);

ci = nlparci(x,resid,'jacobian',J);
nchi2w=(resid).^2./y;
nchi2=sum(nchi2w)./c;
chi2w=(max(yo).*resid).^2./yo;
chi2=sum(chi2w)./c;
resultt=[' t1 = ',num2str(x(2)) ' + ' num2str(x(2)-ci(2,1))...
    ' t2 = ',num2str(x(3)) ' + ' num2str(x(3)-ci(3,1))];
resulta = ['a1 = ',num2str(x(1)) ' + ' num2str(x(1)-ci(1,1)) ' a2 = ',num2str(1-x(1))];
fitresult(1,1)=x(2);
fitresult(1,2)=x(2)-ci(2,1);
fitresult(2,1)=x(3);
fitresult(2,2)=x(3)-ci(3,1);
fitresult(1,3)=x(1);
fitresult(1,4)=x(1)-ci(1,1);
fitresult(2,3)=1-x(1);
fitresult(2,4)=1-x(1);
fitresult(3,1)=x(2)*x(1)+x(3)*(1-x(1));
end 
% fit the result using triple expoential decay
if expn==3

expfun = @(tau,t) real(ifft((fft(irs)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(3)/stepsize)))./tau(3)+tau(2).*exp(-(t/(tau(4)/stepsize)))./tau(4)+(1-tau(1)-tau(2)).*exp(-(t/(tau(5)/stepsize)))./tau(5))))+tau(6)...
    ; %real(ifft((fft(irs)*ones(1,size(y,2))).*

fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));

start_point =[0.5,0.2,16,16,89,background]; 
start_point1 =[0,0,2,5,10,0];
start_point2 =[0.9,0.9,20,50,1000,100];
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',10000,'FunctionTolerance',1e-9);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point,options);
ci = nlparci(x,resid,'jacobian',J);
chi2w=((expfun(x,t)/max(expfun(x,t)))-y).^2./y;
chi2=max(yo).*sum(chi2w)./c;
fitresult(1,1)=x(1)*x(3)+x(2)*x(4)+(1-x(1)-x(2))*x(5)
end
if expn==4
expfun = @(tau,t) real(ifft((fft(irs)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(4)/stepsize)))./tau(4)+tau(2).*exp(-(t/(tau(5)/stepsize)))./tau(5)+tau(3).*exp(-(t/(tau(6)/stepsize)))./tau(6)+(1-tau(1)-tau(2)-tau(3)).*exp(-(t/(tau(7)/stepsize)))./tau(7))))+tau(8)...
    ; %real(ifft((fft(irs)*ones(1,size(y,2))).*
fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));

start_point =[0.1,0.5,0.1,1,10,50,220,0.0001]; 
start_point1 =[0,0,0,0,0,10,20,-10]; 
start_point2 =[0.9,0.4,0.5,10,10,500,1000,100]; 
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',10000,'FunctionTolerance',1e-9);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point,options);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,[],[],options);
ci = nlparci(x,resid,'jacobian',J);
chi2w=((expfun(x,t)/max(expfun(x,t)))-y).^2./y;
chi2=max(yo).*sum(chi2w)./c;

fitresult(1,1)=x(1)*x(4)+x(2)*x(5)+x(3)*x(6)+(1-x(1)-x(2)-x(3))*x(7)
end 

if expn==5
expfun = @(tau,t) real(ifft((fft(irs)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(5)/stepsize)))+tau(2).*exp(-(t/(tau(6)/stepsize)))+tau(3).*exp(-(t/(tau(7)/stepsize)))+tau(4).*exp(-(t/(tau(8)/stepsize)))+(1-tau(1)-tau(2)-tau(3)-tau(4)).*exp(-(t/(tau(9)/stepsize))))))+tau(10)...
    ; %real(ifft((fft(irs)*ones(1,size(y,2))).*
fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));

start_point =[0.5,0.2,0.05,0.05,1,1,5,30,180,background]; 
start_point1 =[0,0,0,0,0,0,0,5,40,-10]; 
start_point2 =[0.9,0.4,0.5,0.5,10,10,10,500,1000,100]; 
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',10000,'FunctionTolerance',1e-9);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point,options);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,[],[],options);
ci = nlparci(x,resid,'jacobian',J);
chi2w=((expfun(x,t)/max(expfun(x,t)))-y).^2./y;
chi2=max(yo).*sum(chi2w)./c;
fitresult(1,2)=1-x(1)-x(2)-x(3)-x(4);
fitresult(1,1)=x(1)*x(5)+x(2)*x(6)+x(3)*x(7)+x(4)*x(8)+(1-x(1)-x(2)-x(3)-x(4))*x(9)
end 

if expn==6
expfun = @(tau,t) real(ifft((fft(irs)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(6)/stepsize)))+tau(2).*exp(-(t/(tau(7)/stepsize)))+tau(3).*exp(-(t/(tau(8)/stepsize)))+tau(4).*exp(-(t/(tau(9)/stepsize)))+tau(5).*exp(-(t/(tau(10)/stepsize)))+(1-tau(1)-tau(2)-tau(3)-tau(4)-tau(5)).*exp(-(t/(tau(11)/stepsize))))))+tau(12)...
    ; %real(ifft((fft(irs)*ones(1,size(y,2))).*
fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));

start_point =[0.5,0.2,0.05,0.02,0.02,1,1,1,1,50,379,background]; 
start_point1 =[0,0,0,0,0,0,0,0,0,5,20,-10]; 
start_point2 =[0.9,0.4,0.5,0.5,0.5,10,10,10,10,500,1000,100]; 
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',10000,'FunctionTolerance',1e-9);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point,options);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,[],[],options);
ci = nlparci(x,resid,'jacobian',J);
chi2w=((expfun(x,t)/max(expfun(x,t)))-y).^2./y;
chi2=max(yo).*sum(chi2w)./c;
fitresult(1,2)=1-x(1)-x(2)-x(3)-x(4)-x(5);
fitresult(1,1)=x(1)*x(6)+x(2)*x(7)+x(3)*x(8)+x(4)*x(9)+x(5)*x(10)+(1-x(1)-x(2)-x(3)-x(4)-x(5))*x(11)

end 

batchfitresult=max(yo)*(expfun(x,t)/max(expfun(x,t)));


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
pos1 = [0.15 0.45 0.8 0.5];
subplot('Position',pos1)
plot(timeaxis,yo);
hold on
plot(timeaxis,max(yo).*expfun(x,t)./max(expfun(x,t)),'LineStyle','-','LineWidth',1.5);
plot(timeaxis,irs./10);
set(gca, 'Fontsize', fs-1)
ylabel('Counts','Fontsize',fs);

 axis([100 400 0.1 75000])
% text(140,1,resultt);
% text(140,0.6,resulta);
legend('Data','Fit','IRF','fontsize',14,'EdgeColor',[1 1 1]);
set(gca, 'YScale', 'log');
pos2 = [0.15 0.14 0.8 0.2];
subplot('Position',pos2)
plot(timeaxis,y-expfun(x,t)./max(expfun(x,t)),'LineStyle','-','LineWidth',1,'Marker','none',...
   'MarkerSize',0.5)
 axis([100 400  -0.2 0.2])
 ylabel('Residuals','Fontsize',fs);

xlabel('Time (ns)','Fontsize',fs);
set(gca,'fontsize',fs-1)
hold off;

figure;
pos1 = [0.15 0.45 0.8 0.5];
subplot('Position',pos1)
plot(timeaxis,y);
hold on
plot(timeaxis,expfun(x,t)./max(expfun(x,t)),'LineStyle','-','LineWidth',1.5);
plot(timeaxis,irs);
set(gca, 'Fontsize', fs-1)
ylabel('Normalized Distribution','Fontsize',fs);

 axis([0 125 0.0006 1.5])
% text(140,1,resultt);
% text(140,0.6,resulta);
legend('Data','Fit','IRF','fontsize',14,'EdgeColor',[1 1 1]);
set(gca, 'YScale', 'log');
pos2 = [0.15 0.14 0.8 0.2];
subplot('Position',pos2)
plot(timeaxis,y-expfun(x,t)./max(expfun(x,t)),'LineStyle','-','LineWidth',1,'Marker','none',...
   'MarkerSize',0.5)
 axis([0 125  -0.2 0.2])
 ylabel('Residuals','Fontsize',fs);

xlabel('Time (ns)','Fontsize',fs);
set(gca,'fontsize',fs-1)
hold off;



