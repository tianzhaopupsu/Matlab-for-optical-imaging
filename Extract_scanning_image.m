close all
%%load file
FIFObfile='C:\Users\Chopin Pro\Desktop\Data\FPGA_confocal_scanning_DATA\Scanning21-10-18_171508';
fid=fopen(FIFObfile,'r','b');
FIFOdata=fread(fid,[3,inf],'int32')';

Scanningimage=zeros(max(FIFOdata(:,1))+1);


for i=1:length(FIFOdata)
    
if (FIFOdata(i,1)>0)&&(FIFOdata(i,2)>0)
Scanningimage(FIFOdata(i,1),FIFOdata(i,2))=Scanningimage(FIFOdata(i,1),FIFOdata(i,2))+FIFOdata(i,3); 
    
end    
end


x = [0 10];
y = [0 10];
imagesc(x,y,Scanningimage);
caxis([10 200]);