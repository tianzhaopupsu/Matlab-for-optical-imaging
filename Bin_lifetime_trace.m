function [binedlf,interlf]=Bin_lifetime_trace(data,binsi)
close all

lfdata=data;
binsize=binsi;
taxis=[0:0.09765625:399.99]';
taxisbin=[0:0.09765625*16:399.99]';
nofp=length(lfdata);
binedlf=zeros(nofp/binsize,1);

for i=binsize:binsize:nofp
  
   for j=1:binsize
   binedlf(i/binsize)=binedlf(i/binsize)+lfdata(i+1-j);
    
   end
end

interlf = interp1(taxisbin,binedlf,taxis);

