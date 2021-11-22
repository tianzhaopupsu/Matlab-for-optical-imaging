close all

FIFOdata=tplfeely10;
reflhis=zeros(4096,1);

photoncounts=length(FIFOdata(:,2));

for i=1:photoncounts
    if FIFOdata(i,2)>0
    reflhis(FIFOdata(i,2))=reflhis(FIFOdata(i,2))+1;
    end
    
end