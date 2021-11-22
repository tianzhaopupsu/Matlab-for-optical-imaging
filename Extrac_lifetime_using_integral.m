close all

dt=aa_timeaxis(2)-aa_timeaxis(1);
maxval=3960;
nrlb=1;
nrhb=130;
a=51;
b=3786;
n=214;
lfres=zeros(n,1);
for i=1:n

ydata=aa_singleinsolusum(a:b,i);

ynoise=ydata-mean(ydata(nrlb:nrhb));
ynorm=ynoise./sum(ynoise);
ar=sum(ynorm);
[M,I] = max(ynorm);
 
y0=ynorm(I:length(ydata));
dy=diff(y0)./dt;
tau=-1./sum(dy);


lfres(i,1)=tau;
end
figure;
plot(y0)


set(gca, 'YScale', 'log')