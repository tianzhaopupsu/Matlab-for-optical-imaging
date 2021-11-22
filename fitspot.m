close all
b=0.2;
fa=importdata('C:\Users\Chopin Pro\Desktop\Scanning_image\NormalizedImages_180821_035145\Z Pixel 0\APD X Left.csv');
x = [0 10];
y = [0 10];
imagesc(x,y,fa);
caxis([0 020]);
colorbar;
c = colorbar;
c.Label.String = 'Counts (kcps)';
set(gca,'fontsize',15);
  xlabel('X (\mum)','FontSize',15)
        ylabel('Y (\mum)','FontSize',15)

map=fa;%load 2D map
map=double(map);
map=medfilt2(map,[3 3]);%median filter, improve the SBR. Make the spots more distinct.
BW=map>80;%value 2 is the threshold of the image
centers = imfindcircles(map,b);%a should be around 0.3
centerround=round(centers);% round the center postion

 sigmaresult=zeros(length(centers),2); %make the storage for the FWHM
for m=1:length(centers)
bi=imcrop(map,[centerround(m,1)-5 centerround(m,2)-5 10 10]);% extract the spot using center postion



  
       amp = max(bi(:));
        xm = length(bi)/2;
        ym = xm;
        sigmax = 0.5;
        sigmay = 0.5;
        offset = 0;
    a0 = [amp xm ym sigmax sigmay offset];  
    opt=optimset('maxfunevals',10000,'maxiter',1e4);
    a = fminsearch(@gaussian2d_fit,a0,opt,bi);
    sigmaresult(m,1)=2*sqrt(2*log(2))*39*(a(4));%find the FWHM for x, 50 is the pixel size
    sigmaresult(m,2)=2*sqrt(2*log(2))*39*(a(5));%find the FWHM for y, 50 is the pixel size
    save('sigmaresult')
    
end
    
%==========================================    
    function chisq = gaussian2d_fit(a,data)
    model = gaussian2d_model(a,data);

    chisq = (model - data).^2;
    chisq = mean(mean(chisq));
    end
%==========================================
function model = gaussian2d_model(a,data)%fit the spot by the 2D Gaussian
    a = abs(a);
    amp = a(1);
    x0 = a(2);
    y0 = a(3);
    sigmax = a(4);
    sigmay = a(5);
    offset = a(6);

    x = (1:length(data(:,1)));
    y = (1:length(data(1,:)));

    for i=1:length(x),
      for j=1:length(y),
        model(i,j) = amp.*exp(-((x(i)-x0).^2)./2./sigmax./sigmax-...
            ((y(j)-y0).^2)./2./sigmay./sigmay )+...
                     offset;%the model of 2D Gaussian

      end
    end


end
