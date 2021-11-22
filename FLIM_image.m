close all

filefolder='\\192.168.1.101\newdata\tian\Y2021\M11\D19\FLIM_191121_164543';
binsi=16;
imagesize=256;

taxis=[0:0.09765625:399.99]';
taxisbin=[0:0.09765625*16:399.99]';
filePattern = fullfile(filefolder, '*');
matFiles = dir(filePattern);
filenum=length(matFiles);
FLIMimag=zeros(imagesize,imagesize);
h = waitbar(0,'Processing FLIM image...');
for i=3:filenum

xpixel=str2double(matFiles(i).name(20:22));
ypixel=str2double(matFiles(i).name(23:25));
if xpixel>0 && xpixel<(imagesize+1)&& ypixel>0 && ypixel<(imagesize+1)
f = fullfile(filefolder,matFiles(i).name);
fid = [importdata(f)]';
[binedlf,interlf]=Bin_lifetime_trace(fid,binsi);
[tau]=lifetime_fit_tri(interlf,taxis);
%[exptaufinal]=Integral_lifetime_func(interlf,taxis);
FLIMimag(xpixel,ypixel)=tau;
end
waitbar(i/filenum,h)
end