function [tau]=lifetime_fit_tri(data,tax)

%cw is the center postion of the instrument respons function
%y is the fluoresence decay of the sample
%how many exponentials you want in this fitting, 1==1exp 2==2exp 3==3exp.
%setpsize is the time interval between pixels.
close all

a=970;
b=3900;
yo=data(a:b);
timeaxis=tax(a:b);

background=0;
%yo=yc(5:3897,:);
%yo=yc(53:3897,:);


stepsize=timeaxis(2)-timeaxis(1);
%cw=63;
cw=310;
%cw=23;
FWHM=0.8;%0.944  0.34 in ns
width=FWHM/stepsize;

c=length(yo);
t=[0:(c-1)]';
y=yo./max(yo);
y(y==0)=1;

irs=exp(-((t-cw)/width).^2);



expfun = @(tau,t) real(ifft((fft(irs)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(3)/stepsize)))./tau(3)+tau(2).*exp(-(t/(tau(4)/stepsize)))./tau(4)+(1-tau(1)-tau(2)).*exp(-(t/(tau(5)/stepsize)))./tau(5))))+tau(6)...
    ; %real(ifft((fft(irs)*ones(1,size(y,2))).*

fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));

start_point =[0.5,0.2,16,16,89,background]; 
start_point1 =[0,0,2,5,50,-10];
start_point2 =[0.9,0.9,20,50,300,100];
options = optimoptions('lsqnonlin','Display','off','MaxFunctionEvaluations',100,'FunctionTolerance',1e-9);
[x] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);

tau=x(1)*x(3)+x(2)*x(4)+(1-x(1)-x(2))*x(5);
end


