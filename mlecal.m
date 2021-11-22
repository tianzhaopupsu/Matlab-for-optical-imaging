function [Dc, sigma, min_var_index, sigma_Dc] = mlecal
close all

%fid = fopen('\\192.168.1.101\newdata\tian\Y2021\M09\D06\Tracking_trajactories _matlab\1');
fid = fopen('C:\Users\Chopin Pro\Desktop\Data\tracking21-11-08_170650');
chan1=fread(fid,[11, inf],'double');
chan1=chan1';
trackDataFilt=zeros(length(chan1),11);
trackDataFilt(:,1)=chan1(:,1);
trackDataFilt(:,2)=chan1(:,2);
trackDataFilt(:,3)=chan1(:,3);
trackDataFilt(:,4)=chan1(:,4);
trackDataFilt(:,5)=chan1(:,5);
trackDataFilt(:,6)=chan1(:,6);
trackDataFilt(:,7)=chan1(:,7);
trackDataFilt(:,8)=chan1(:,8);
trackDataFilt(:,9)=chan1(:,9);
trackDataFilt(:,10)=chan1(:,10);
trackDataFilt(:,11)=chan1(:,11);
a=length(trackDataFilt);
trackDataFilt=trackDataFilt(0.1*a:0.90*a,1:11);
intensity=trackDataFilt(:,7)+trackDataFilt(:,8)+trackDataFilt(:,9)+trackDataFilt(:,10)+trackDataFilt(:,11);
xyz=[trackDataFilt(:,1) trackDataFilt(:,2) trackDataFilt(:,3)];
xyz=xyz;
n=length(xyz);
sample_dt = 1e-3; % coordinate data logging time
bin_time = 10e-6; % bin time per photon-counting sample used in Kalman filter
lambda = 0.003; % the AR parameter used in Kalman filter
% eta = 0.0068559; %Ns/m2
% T = 298.15; % 20 C
% kB = 1.38e-23; % the Boltzmann constant
% D_trans_ideal = kB.*T./(6.*pi.*eta.*100e-9); % theoretical translation diff.


% Set up anonymous function for the Kalman filter correction factor [1].
A = @(m) (2-2.*(1-lambda).^(1+m)-2.*lambda-2.*m.*lambda+m.*lambda.^2)./...
         (m.*(-2+lambda).*lambda);

% Set up anonymous function for fitting noise model [1,3].
Dcsq = @(parm,D_app,dta) mean((D_app-parm(1)-parm(2)./dta).^2);

numberOfDeltaT=floor(n/4);
N=zeros(numberOfDeltaT,1);
D_app = zeros(numberOfDeltaT,3); 
tic
for dt = 1:numberOfDeltaT
   deltaCoords = xyz(1+dt:end,:) - xyz(1:end-dt,:);
   squaredDisplacement = deltaCoords.^2; %# dx^2+dy^2+dz^2
   N(dt)=size(squaredDisplacement,1);
   D_app(dt,:) = mean(squaredDisplacement)./dt./sample_dt./2;   
   %waitbar(dt/numberOfDeltaT,hWait)
end

delta = (1:numberOfDeltaT)'.*sample_dt;
m = delta / bin_time; % number of photon-counting samples within this time lag

Dx_app=D_app(:,1);
Dy_app=D_app(:,2);
Dz_app=D_app(:,3);

% fit the scaled D_app to a noise model for the corrected diff. coef., Dc
if length(delta) < 40
    L=length(delta);
else
    L=length(delta);
end


dta = delta(1:L);
% -- calculate x
D_app_tmp = Dx_app(1:L)./A(m(1:L));
parm_opt = fminsearch( @(parm) Dcsq(parm,D_app_tmp,dta),[D_app_tmp(end), 0]);
% Dcx = parm_opt(1);
sigmax = sqrt(parm_opt(2)); % spatial localization prevision of the x-axis
var_Dcx = 2.*(parm_opt(2) + delta.*parm_opt(1)).^2./N./lambda./delta.^2;

% -- calculate y
D_app_tmp = Dy_app(1:L)./A(m(1:L));
parm_opt = fminsearch( @(parm) Dcsq(parm,D_app_tmp,dta),[D_app_tmp(end), 0]);
% Dcy = parm_opt(1);
sigmay = sqrt(parm_opt(2)); % spatial localization prevision of the y-axis
var_Dcy = 2.*(parm_opt(2) + delta.*parm_opt(1)).^2./N./lambda./delta.^2;

% -- calculate z
D_app_tmp = Dz_app(1:L)./A(m(1:L));
parm_opt = fminsearch( @(parm) Dcsq(parm,D_app_tmp,dta),[D_app_tmp(end), 0]);
% Dcz = parm_opt(1);
sigmaz = sqrt(parm_opt(2)); % spatial localization prevision of the z-axis
var_Dcz = 2.*(parm_opt(2) + delta.*parm_opt(1)).^2./N./lambda./delta.^2;

% -- find the optimal lag time for minimal D uncertainty
[tmp,min_varx_index] = min(var_Dcx);
[tmp,min_vary_index] = min(var_Dcy);
[tmp,min_varz_index] = min(var_Dcz);
% -- calculate the uncertainty at the minimal uncertainty time lag
sigma_Dcx = sqrt(var_Dcx(min_varx_index));
sigma_Dcy = sqrt(var_Dcy(min_vary_index));
sigma_Dcz = sqrt(var_Dcz(min_varz_index));
% -- calculate the diffusion coefficient @ minimal uncertainty time lag
Dcx = Dx_app(min_varx_index)./A(m(min_varx_index)) - ...
      sigmax^2/delta(min_varx_index);
Dcy = Dy_app(min_vary_index)./A(m(min_vary_index)) - ...
      sigmay^2/delta(min_vary_index);
Dcz = Dz_app(min_varz_index)./A(m(min_varz_index)) - ...
      sigmaz^2/delta(min_varz_index);
  
  
Dc=[Dcx,Dcy,Dcz];
sigma=[sigmax,sigmay,sigmaz];
min_var_index=[min_varx_index,min_vary_index,min_varz_index];
sigma_Dc=[sigma_Dcx,sigma_Dcy,sigma_Dcz];
toc


delta_finer = (delta(1):sample_dt/5:delta(end))';
Rz(:,1)=Dx_app;
Rz(:,2)=Dy_app;
Rz(:,3)=Dz_app;
Rz(:,4)=Dx_app./A(m);
Rz(:,5)=Dy_app./A(m);
Rz(:,6)=Dz_app./A(m);
Rz(:,7)=delta*1e3;
Rz1(:,1)=Dcx+sigmax^2./delta_finer;
Rz1(:,2)=Dcy+sigmay^2./delta_finer;
Rz1(:,3)=Dcz+sigmaz^2./delta_finer;
Rz1(:,4)=delta_finer*1e3;
x=xyz(:,1);
y=xyz(:,2);
z=xyz(:,3);
DT=(Dcx+Dcy+Dcz)/3;
sigmaDT=sqrt(sigma_Dcx.^2+sigma_Dcy.^2+sigma_Dcz.^2)/3;
disp([num2str(DT) ',' num2str(sigmaDT)]);
x=x-x(1,1);
y=y-y(1,1);
z=z-z(1,1);
%%intensity=(intensity/0.4);
 x=x;
c=1:length(trackDataFilt);
c=c./1000;
fs=15;
h = surface(...
  'XData',[x(:) x(:)],...
  'YData',[y(:) y(:)],...
  'ZData',[z(:) z(:)],...
  'CData',[c(:) c(:)],...
  'FaceColor','none',...
  'EdgeColor','flat',...
  'Marker','none');
disp(['std: x:' num2str(std(x)) ', y:' num2str(std(y)) ', z:' num2str(std(z)) ', time:' num2str(a/1000)] )
set(gca,'fontsize',fs-1)
     xlabel('X (\mum)','FontSize',fs)
        ylabel('Y (\mum)','FontSize',fs)
        zlabel('Z (\mum)','FontSize',fs)
        c = colorbar;
    
        
        c.Location = 'southoutside';
       set(c,'fontsize',fs-1);
c.Label.String = ['Time (s)'];
       %c.Label.String = ['Trajectory length = ' num2str(round(n/1000,1)) 's'];%'Trajectory length = 23.5s'
 grid on
set(gca,'GridLineStyle','--')        
figure;
        plot(intensity);
        set(gca,'fontsize',15)
        xlabel('Time (ms)','FontSize',15)
        ylabel('Total Photon Counts (kcps)','FontSize',15)


figure
subplot(3,1,1);
plot(delta*1e3, Dx_app./A(m), 'ko', ...
     delta_finer*1e3, Dcx+sigmax^2./delta_finer, 'r-');
yl = ylim;
axis([1 140 yl(1) yl(2)]);
txt = ['D_x = ' sprintf('%.3f',Dcx) ' \pm ' sprintf('%.3f',sigma_Dcx) ...
       ' \mum^2/s'];
text(1, sum(yl)/2, txt);
legend( 'x' );
ylabel( 'D_{app} (\mum^2/s)', 'Fontsize',8);
%title('Apparent diffusion coeff.' );%trans\_rot\_temperature.m: Fig. 3, apparent diffusion coeff.

subplot(3,1,2);
plot(delta*1e3, Dy_app./A(m), 'ko', ...
     delta_finer*1e3, Dcy+sigmay^2./delta_finer, 'r-');
yl = ylim;
axis([1 140 yl(1) yl(2)]);
txt = ['D_y = ' sprintf('%.3f',Dcy) ' \pm ' sprintf('%.3f',sigma_Dcy) ...
       ' \mum^2/s'];
text(1, sum(yl)/2, txt);
set(gca,'Fontsize',7)
legend( 'Experimental Data' ,'Fit', 'Fontsize',7);
ylabel( 'D_{app} (\mum^2/s)' , 'Fontsize',8);

subplot(3,1,3);
plot(delta*1e3, Dz_app./A(m), 'ko', ...
     delta_finer*1e3, Dcz+sigmaz^2./delta_finer, 'r-');
yl = ylim;
axis([1 140 yl(1) yl(2)]);
txt = ['D_z = ' sprintf('%.3f',Dcz) ' \pm ' sprintf('%.3f',sigma_Dcz) ...
       ' \mum^2/s'];
text(1, sum(yl)/2, txt);
legend( 'z' );
ylabel( 'D_{app} (\mum^2/s)' );
xlabel( 'lag time (\delta, ms)' );
%%subplot(4,1,4);
figure;
plot(delta*1e3,Dx_app,'r',delta*1e3, Dy_app,'b', delta*1e3,Dz_app,'g');
yl = ylim;
% axis([0 500 0 10]);
xlabel( 'time lag (ms)' );
ylabel( 'D_{app} (\mum^2/s)' );
legend('x','y','z')
title('trans\_rot\_temperature.m: Fig. 2, evaluate consistency' );
saveas(gca,'D_c.fig');
save('Dcpamrs','Dc', 'sigma', 'min_var_index', 'sigma_Dc','Rz','Rz1','sigmax','sigmay','sigmaz');

