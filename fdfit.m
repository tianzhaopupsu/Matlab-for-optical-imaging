
%cw is the center postion of the instrument respons function
%y is the fluoresence decay of the sample
%how many exponentials you want in this fitting, 1==1exp 2==2exp 3==3exp.
%setpsize is the time interval between pixels.
close all
yc=(aas35_500ns(3:3997,3)./aas35_500ns(3800,3))-(aas35_500ns(3:3997,6)./aas35_500ns(3800,6));
expn=1;
irff=2;
background=0.000001;
%yo=yc(5:3897,:);
%yo=yc(53:3897,:);
yo=yc+2;

funda=aafund(5:3897,2);
irf=funda;
stepsize=0.122;
%cw=63;
cw=71;
FWHM=1.2;%0.944  0.34 in ns
width=FWHM/stepsize;
fs=8; %font size
c=length(yo);
t=[0:(c-1)]';
y=yo./max(yo);
fitresult=zeros(4,4);

if irff==1
irf=irf./max(irf);
    irs=irf;

end

if irff==2
irs=exp(-((t-cw)/width).^2);
end

% fit the result using single expoential decay
if expn==1

expfun = @(tau,t) (real(ifft((fft(irs)*ones(1,size(y,2))).*fft(exp(-(t/(tau(1)/stepsize)))+tau(2)))));

fitres =  @(tau) sum((((expfun(tau,t)/max(expfun(tau,t)))-y).^2)./y);
fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y)).^2./y);

start_point =[43 background*160000]; 
start_point1 =[12 background*0]; 
start_point2 =[102 background*150000]; 
B = fminsearch(fitres, start_point);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2);
ci = nlparci(x,resid,'jacobian',J);
nchi2w=(resid).^2./y;
nchi2=sum(nchi2w)./c;
chi2w=(max(yo).*resid).^2./yo;
chi2=sum(chi2w)./c;
resultt=['t = ',num2str(x(1)) '+' num2str(x(1)-ci(1,1))];
resulta = num2str(1);
fitresult(1,1)=x(1);
fitresult(1,2)=x(1)-ci(1,1);
end 

% fit the result using double expoential decay
if expn==2

expfun = @(a,t) (real(ifft((fft(irs)*ones(1,size(y,2))).*fft(a(1)*exp(-(t/(a(2)/stepsize)))+(1-a(1))*exp(-(t/(a(3)/stepsize)))+a(4)))));

fitres =  @(tau) sum((((expfun(tau,t)/max(expfun(tau,t)))-y).^2)./y);
fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y)).^2./y);
start_point =[0.4,35,25,background*12]; 
B = fminsearch(fitres, start_point);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,[],[]);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point./2,start_point.*2);

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

expfun = @(tau,t) real(ifft((fft(irs)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(2)/stepsize)))+tau(3).*exp(-(t/(tau(4)/stepsize)))+tau(5).*exp(-(t/(tau(6)/stepsize)))+tau(7))))...
    ; %real(ifft((fft(irs)*ones(1,size(y,2))).*

fitres =  @(tau) sum((((expfun(tau,t)/max(expfun(tau,t)))-y).^2)./y);
fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y)).^2./y);

start_point =[0.058,20,0.041,280,0.02,600,background]; 
start_point1 =[0.048,10,0.021,10,0.006,200,background]; 
start_point2 =[0.068,50,0.051,330,0.03,1310,background]; 
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',10000,'FunctionTolerance',1e-9);
B = fminsearch(fitres, start_point);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);
%[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point,start_point,options);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,[],[],options);
ci = nlparci(x,resid,'jacobian',J);
nchi2w=(resid).^2./y;
nchi2=sum(nchi2w)./c;
chi2w=(max(yo).*resid).^2./yo;
chi2=sum(chi2w)./c;
resultt=['t1 = ',num2str(x(2)) '+' num2str(x(2)-ci(2,1)) ' t2 = ',num2str(x(4)) '+' num2str(x(4)-ci(4,1))...
    ' t3 = ',num2str(x(6)) '+' num2str(x(6)-ci(6,1))];
sumpor = x(1)+x(3)+x(5);
resulta = ['t1 = ',num2str(x(1)/sumpor) '+' num2str((x(1)-ci(1,1))/sumpor) ' t2 = ',num2str(x(3)/sumpor) '+' num2str((x(3)-ci(3,1))/sumpor)...
    ' t3 = ',num2str(x(5)/sumpor) '+' num2str((x(5)-ci(5,1))/sumpor)];
fitresult(1,1)=x(2);
fitresult(1,2)=x(2)-ci(2,1);
fitresult(2,1)=x(4);
fitresult(2,2)=x(4)-ci(4,1);
fitresult(3,1)=x(6);
fitresult(3,2)=x(6)-ci(6,1);
fitresult(1,3)=x(1)/sumpor;
fitresult(1,4)=(x(1)-ci(1,1))/sumpor;
fitresult(2,3)=(x(3))/sumpor;
fitresult(2,4)=(x(3)-ci(3,1))/sumpor;
fitresult(3,3)=(x(5))/sumpor;
fitresult(3,4)=(x(5)-ci(5,1))/sumpor;
fitresult(4,1)=x(2)*x(1)/sumpor+x(4)*x(3)/sumpor+x(6)*x(5)/sumpor;


end 

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
plot(y);
hold on
plot(expfun(x,t)/max(expfun(x,t)));
set(gca, 'Fontsize', fs-1)
ylabel('Normalized Distribution','Fontsize',fs);

xlabel('Time (ns)','Fontsize',fs);
 axis([0 3990 0.001 1.5])
text(140,1.2,resultt);
text(140,1,resulta);

set(gca, 'YScale', 'log');
hold off;
expfunr=expfun(x,t);
aexpfunrn=expfun(x,t)/max(expfun(x,t));




