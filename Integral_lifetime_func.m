function [exptaufinal]=Integral_lifetime_func(data,ta)
close all

expdata=data;
tax=ta;
nrlb=1000;
nrhb=1100;
a=970;
b=3900;
numofw=1;
expdata(expdata==0)=1;
for i=1:numofw

%ydata=aa_singleinsolusum(a:b,i);

ydata=expdata(a:b,i);

[M,I] = max(ydata);
ynoise=ydata-mean(ydata(nrlb:nrhb));
y0=ynoise(I:length(ydata));
y0( y0 <= 0 ) = 1;
tp=tax(I:length(ydata));
exptau=0;
for m=2:length(y0)
taucontri=(tp(m)-tp(1))*y0(m);   

exptau=exptau+taucontri;

end
exptaufinal=exptau/sum(y0);

end
