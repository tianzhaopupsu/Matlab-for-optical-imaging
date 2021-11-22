f=0.18;  %focal length of the lens in mm
Power=0.013; %micro W
rep=0.2*10^6;
wavelength=400*10^-7;
r=f*wavelength/(pi*3); % 3 in mm
areaobj=pi*r.^2; %cm^2
pulsenergy=0.83.*Power/rep;
density=pulsenergy./areaobj;