close all
t=[0:0.030517578125:124.99]';
dt=t(2)-t(1);
a1=1;
a2=0;
a3=0;
a4=0;
a5=0;
t1=3.85;
t2=1;
t3=1;
t4=1;
t5=1;
mt=2;
n=500;
maxval=length(t);
fs=15;
irf=exp(-((t-20)/0.1).^2);
irf=irf/sum(irf);
if mt==1
ymod=(a1*exp(-t/t1)/t1+a2*exp(-t/t2)/t2+a3*exp(-t/t3)/t3+a4*exp(-t/t4)/t4+a5*exp(-t/t5)/t5);
end
if mt==2
ymod=(a1*exp(-t/t1)+a2*exp(-t/t2)+a3*exp(-t/t3)+a4*exp(-t/t4)+a5*exp(-t/t5));
end
ycolv=conv( ymod, irf);
ycolv=ycolv(1:maxval);
ycolv=ycolv./sum(ycolv);
ltsresult=zeros(n,1);
for i= 1:n

ynoise=poissrnd(ycolv*1000000+100);
ynoise=ynoise-mean(ynoise(3563:4050));

%ynorm=ycolv./sum(ycolv);
%ynorm=ymod./sum(ymod);
ynorm=ynoise./sum(ynoise);
b=sum(ynorm);


 [M,I] = max(ynorm);
 [IM,II] = max(irf);
 
 y0=ynorm(I:maxval);
% irfarea=irf(II:maxval);
% irfareanorm=irfarea./sum(irfarea);
% dirf=diff(irfareanorm)./diff(t(II:maxval));
% tauirf=-1./sum(dirf);
% 
% y0=fft(ycolv-100)./fft(irf);
% y0=y0*ones(1,size(ycolv,2));
% y1=real(ifft(y0));
% ynorm=y1./sum(y1);
dy=diff(y0)./dt;
tau=-1./sum(dy);
ddy=diff(dy)./dt;
tau1=1/sum(ddy);
% sumtau=zeros(2*I);
% for i=1:2*I
% spoint=i;
% y0=ycolv(spoint:100001)-100;
% ynorm=y0./sum(y0);
% dy=diff(ynorm)./diff(t(spoint:100001));
% tau=-1./sum(dy);
% sumtau(i)=tau-tauirf;
% 
% end
abs_tau=a1*t1+a2*t2+a3*t3+a4*t4+a5*t5;
abs_sum=a1/t1+a2/t2+a3/t3+a4/t4+a5/t5;
abs_tau1=t1*a1/t1/abs_sum+t2*a2/t2/abs_sum+t3*a3/t3/abs_sum+t4*a4/t4/abs_sum+t5*a5/t5/abs_sum;
abs_tau2=t1*t2/(a1*t2+a2*t1);

ltsresult(i,1)=tau;
end
avg=mean(ltsresult)
dev=std(ltsresult)
plot(t,ynoise)
hold on

%plot(t,irf)

set(gca, 'YScale', 'log')
set(gca,'fontsize',fs-1)
ylabel('Counts','Fontsize',fs);
xlabel('Time (ns)','Fontsize',fs);
%axis([0 6000  1 100])