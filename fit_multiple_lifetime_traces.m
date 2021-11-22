
close all
a=61;
b=3786;
c=0;
datamatrix=aa_interface_airfix;
funda=aa_irf((a-c:b-c),1);

timeaxis=aa_timeaxis(a:b);


%%%%--------tetra initial
% a0 = [0.5, 0.2, 0.1];            
% tau0 = [1, 1, 50, 400];
% bkg = 0;
% start_point = [a0, tau0, bkg];
% start_point1 =[0,0,0,0,0,10,50,-10]; 
% start_point2 =[0.9,0.4,0.5,10,10,500,1000,100]; 


%%%%--------penta initial


a0 = [0.5, 0.2, 0.05,0.05];            
tau0 = [1, 1, 5, 10, 80];
bkg = 0;
start_point = [a0, tau0, bkg];
start_point1 =[0,0,0,0,0,0.5,6,5,50,-10]; 
start_point2 =[0.9,0.6,0.6,0.6,5,80,80,100,1000,100]; 

stepsize=0.30517578125;
[m,n]=size(datamatrix);
batchresult=zeros(n,12);
batchfitresult=zeros(b-a+1,n);
irf=funda/sum(funda);
t=[0:(length(funda)-1)]';
for i = 1:n

yo=datamatrix((a:b),i);
c=length(yo);
y=yo./max(yo);
%%%%%%%%%tetra-exp
% 
% expfun = @(tau,t) real(ifft((fft(irf)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(4)/stepsize)))+tau(2).*exp(-(t/(tau(5)/stepsize)))+tau(3).*exp(-(t/(tau(6)/stepsize)))+(1-tau(1)-tau(2)-tau(3)).*exp(-(t/(tau(7)/stepsize))))))+tau(8)...
%     ; %real(ifft((fft(irs)*ones(1,size(y,2))).*
% fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));
% options = optimoptions('lsqnonlin','MaxFunctionEvaluations',10000,'FunctionTolerance',1e-9);
% [x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);
% ci = nlparci(x,resid,'jacobian',J);
% chi2w=((expfun(x,t)/max(expfun(x,t)))-y).^2./y;
% chi2=max(yo).*sum(chi2w)./c;
% batchresult(i,1:8)=x;
% batchresult(i,9)=chi2;
% batchresult(i,10)=x(1)*x(4)+x(2)*x(5)+x(3)*x(6)+(1-x(1)-x(2)-x(3))*x(7);
% batchfitresult(:,i)=max(yo)*(expfun(x,t)/max(expfun(x,t)));


%%%%%%%%%penta-exp

expfun = @(tau,t) real(ifft((fft(irf)*ones(1,size(y,2))).*fft(tau(1).*exp(-(t/(tau(5)/stepsize)))+tau(2).*exp(-(t/(tau(6)/stepsize)))+tau(3).*exp(-(t/(tau(7)/stepsize)))+tau(4).*exp(-(t/(tau(8)/stepsize)))+(1-tau(1)-tau(2)-tau(3)-tau(4)).*exp(-(t/(tau(9)/stepsize))))))+tau(10)...
    ; 
fitresls =  @(tau) ((((expfun(tau,t)/max(expfun(tau,t)))-y))./sqrt(y));
options = optimoptions('lsqnonlin','MaxFunctionEvaluations',10000,'FunctionTolerance',1e-9);
[x,resnorm,resid,exitflag,output,lambda,J] = lsqnonlin(fitresls, start_point,start_point1,start_point2,options);
ci = nlparci(x,resid,'jacobian',J);
chi2w=((expfun(x,t)/max(expfun(x,t)))-y).^2./y;
chi2=max(yo).*sum(chi2w)./c;
batchresult(i,1:10)=x;
batchresult(i,11)=chi2;
batchresult(i,12)=x(1)*x(5)+x(2)*x(6)+x(3)*x(7)+x(4)*x(8)+(1-x(1)-x(2)-x(3)-x(4))*x(9);
batchresult(i,13)=(1-x(1)-x(2)-x(3)-x(4));
batchfitresult(:,i)=max(yo)*(expfun(x,t)/max(expfun(x,t)));



fprintf( 1, ['Finished the analysis for the lifetime trace No.' num2str(i)] );
end 


figure;
plot(batchresult(:,11))

