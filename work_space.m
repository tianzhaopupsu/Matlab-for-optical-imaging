                                                                                n=23;
spec2=importdata('C:\Users\Chopin Pro\Desktop\Lifetime\D07282021_bg.asc');
% spec3=importdata('C:\Users\Chopin Pro\Desktop\Lifetime\D07282021_irf500nsa.asc');
% spec4=importdata('C:\Users\Chopin Pro\Desktop\Lifetime\D07282021_irf50nsa.asc');
spec1(:,n*3+1)=spec2.data;
spec1(:,n*3+2)=spec3.data;
spec1(:,n*3+3)=spec4.data;