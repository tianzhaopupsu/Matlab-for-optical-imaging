
a=0.055270;
b=10;
% eta = 0.0010049; %Ns/m2 water
% eta = 0.0083650; %Ns/m2 1:1 w/g
%eta = 0.0103650; %Ns/m2 1:1.15 w/g
%eta = 0.016062; %Ns/m2 1:1.5 w/g
%eta = 0.025249;  %Ns/m2 1:2 w/g
D_trans_ideal = 1.38*10^-23*293.15/(6.*pi.*a.*b*10^-9)/10^-12


radius=1.38*10^-23*293.15/(D_trans_ideal*10^-12*6.*pi.*a.*10^-9);



aa_diameter0916=1*1.38*10^-23*293.15./(aa_mleresult0916(:,1)*10^-12*6.*pi.*a.*10^-9);