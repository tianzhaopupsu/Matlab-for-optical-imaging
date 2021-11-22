function [dsresult] = batch_mle
close all
%filenum=52;



filePattern = fullfile('\\192.168.1.101\newdata\tian\Y2021\M09\D16\Tracking_09162021', '*');
matFiles = dir(filePattern);
filenum=length(matFiles);
dsresult=zeros(filenum,3);
for i = 3:filenum

f = fullfile('\\192.168.1.101\newdata\tian\Y2021\M09\D16\Tracking_09162021',matFiles(i).name);
fid = fopen(f);
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


DT=(Dcx+Dcy+Dcz)/3;
sigmaDT=sqrt(sigma_Dcx.^2+sigma_Dcy.^2+sigma_Dcz.^2)/3;
dsresult(i,1)=DT;
dsresult(i,2)=sigmaDT;
dsresult(i,3)=a/1000;
end

