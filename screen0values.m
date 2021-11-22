close all

aa_fix=aa_800;
[m,n]=size(aa_fix);
for i = 1:n

    for j=61:3786
        
        if aa_fix(j,i)<1
            
            aa_fix(j,i)= aa_fix(j-1,i);
        
        end
        
    end
    
    
   
    
end