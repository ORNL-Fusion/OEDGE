for n=1:16
    name=strcat('chi_',num2str(n,'%.3d'));
    name=strcat(name,'.txt');
    chi=load(name);
    [m,p]=size(chi);
    for i=1:m
        l=0;
        waitingforchi1=1;
        for j=1:p
           if((chi(i,j)==1)&&(waitingforchi1==1))
               if(l==1)
                  chi(i,j-1)=1; 
               end
               if(l==2)
                    chi(i,j-1)=1;
                    chi(i,j-2)=1;
               end
               waitingforchi1=0;
               l=0;
           end
           if(chi(i,j)==0)
               l=l+1;
               waitingforchi1=1;
           end
        end
        if((chi(i,p)==0)&&(chi(i,p-1)==1))
            chi(i,p)=1;
        end
    end
    save(name,'chi','-ascii');
end