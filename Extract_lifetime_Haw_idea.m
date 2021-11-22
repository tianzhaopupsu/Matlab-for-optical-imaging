close all

expdata=interlf;
tax=taxis;
nrlb=1000;
nrhb=1100;
a=970;
b=3900;
numofw=1;
expdata(expdata==0)=1;
lfres=zeros(numofw,1);
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
lfres(i,1)=exptaufinal;
fprintf( 1, [' Traj. No. 1 Lifetime = ' num2str(exptaufinal) ' ns;'] );
end
figure;
plot(y0)


set(gca, 'YScale', 'log')