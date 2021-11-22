close all
t=[0:0.030517578125:124.99]';
dt=0.030517578125;
a1=0.073;
a2=0.367;
a3=0.372;
a4=0.188;
a5=0;
t1=0.0345;
t2=4.25;
t3=22.586;
t4=0.633;
t5=1;
mt=1;
n=500;
maxval=length(t);
fs=15;

if mt==1
ymod=(a1*exp(-t/t1)/t1+a2*exp(-t/t2)/t2+a3*exp(-t/t3)/t3+a4*exp(-t/t4)/t4+a5*exp(-t/t5)/t5);
end
if mt==2
ymod=(a1*exp(-t/t1)+a2*exp(-t/t2)+a3*exp(-t/t3)+a4*exp(-t/t4)+a5*exp(-t/t5));
end
abs_tau=a1*t1+a2*t2+a3*t3+a4*t4+a5*t5;
abs_sum=a1/t1+a2/t2+a3/t3+a4/t4+a5/t5;
abs_tau1=t1*a1/t1/abs_sum+t2*a2/t2/abs_sum+t3*a3/t3/abs_sum+t4*a4/t4/abs_sum+t5*a5/t5/abs_sum;
resultss=zeros(57,2);
for wn=10:10
width=wn/23.55;

irf=exp(-((t-20)/(1.44*width)).^2);
irf=irf/sum(irf);
ycolv=conv( ymod, irf);
ycolv=ycolv(1:maxval);
ltsresult=zeros(10,1);
for sn=1:1
ynoise=poissrnd(ycolv*1000000+100);
[M,I] = max(ynoise);
ynoise=ynoise-mean(ynoise(3563:4050));
y0=ynoise(I:maxval);
y0( y0 <= 0 ) = 0;
tp=t(I:maxval);
exptau=0;
for m=2:length(y0)

 taucontri=(tp(m)-tp(1))*y0(m);   


exptau=exptau+taucontri;

end
exptaufinal=exptau/sum(y0);
ltsresult(sn,1)=exptaufinal;
end

resultss(wn-3,1)=mean(ltsresult);
resultss(wn-3,2)=std(ltsresult);

end
FWHMMaxis=[4:1:60];
%plot(FWHMMaxis/23.55,resultss(:,1))


plot(t,ynoise)

set(gca, 'YScale', 'log')
set(gca,'fontsize',fs-1)
ylabel('Counts','Fontsize',fs);
xlabel('Time (ns)','Fontsize',fs);

figure;
errorbar(FWHMMaxis/10,resultss(:,1),resultss(:,2))

set(gca,'fontsize',fs-1)
ylabel('Lifetime (ns)','Fontsize',fs);
xlabel('FWHM (ns)','Fontsize',fs);

