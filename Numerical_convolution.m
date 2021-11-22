close all
t=[0:0.030517578125:124.99]';
dt=0.030517578125;
a1=1;
a2=0;
a3=0;
a4=0;
a5=0;

t2=1;
t3=1;
t4=1;
t5=1;
mt=1;
n=200;
maxval=length(t);
fs=15;
sre1=zeros(n,n);
sre2=zeros(n,n);
sre3=zeros(n,n);
sre4=zeros(n,n);

for i=1:n
    
t1=i/10;
for j=1:n
FWHM=j*0.015;
sigma=FWHM/2.355;
sigma=sigma*1.44;
irf=exp(-((t-20)/sigma).^2);
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



[M,I] = max(ycolv);

y0=ycolv(I:maxval);

tp=t(I:maxval);
exptau=0;
for m=2:length(y0)

 taucontri=(tp(m)-tp(1))*y0(m);   


exptau=exptau+taucontri;

end
exptaufinal=exptau/sum(y0);
sre1(i,j)=t1;
sre2(i,j)=FWHM;
sre3(i,j)=exptaufinal;
sre4(i,j)=(1-(absb(exptaufinal-t1)/t1))*100;
end
i
end
plot(t,ynorm)
hold on

%plot(t,irf)

set(gca, 'YScale', 'log')
set(gca,'fontsize',fs-1)
ylabel('Counts','Fontsize',fs);
xlabel('Time (ns)','Fontsize',fs);
%axis([0 6000  1 100])