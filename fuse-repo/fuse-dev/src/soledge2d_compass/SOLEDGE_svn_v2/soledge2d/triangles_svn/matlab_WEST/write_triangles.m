fid=fopen('vert2knots.txt','w');
fprintf(fid,strcat(num2str(ntri),'\t\t!number of knots\n'));
for i=1:ntri
    knot(i).vertnum=0;
    knot(i).inum=[];
    knot(i).jnum=[];
    knot(i).knum=[];
    for j=1:16
        Nx=zone(j).Nx;
        Nz=zone(j).Nz;
        
        [m,n]=find(zone(j).ntrinum==i);
        if(not(isempty([m,n])))
            knot(i).R=zone(j).Rtri(m,n);
            knot(i).Z=zone(j).Ztri(m,n);
        end
        if(m-1>=1)
            if(n-1>=1)
                knot(i).vertnum=knot(i).vertnum+1;
                knot(i).inum=[knot(i).inum,m-1];
                knot(i).jnum=[knot(i).jnum,n-1];
                knot(i).knum=[knot(i).knum,j];
            end
            if(n<=Nz)
                knot(i).vertnum=knot(i).vertnum+1;
                knot(i).inum=[knot(i).inum,m-1];
                knot(i).jnum=[knot(i).jnum,n];
                knot(i).knum=[knot(i).knum,j];
            end
        end
        if(m<=Nx)
            if(n-1>=1)
                knot(i).vertnum=knot(i).vertnum+1;
                knot(i).inum=[knot(i).inum,m];
                knot(i).jnum=[knot(i).jnum,n-1];
                knot(i).knum=[knot(i).knum,j];
            end
            if(n<=Nz)
                knot(i).vertnum=knot(i).vertnum+1;
                knot(i).inum=[knot(i).inum,m];
                knot(i).jnum=[knot(i).jnum,n];
                knot(i).knum=[knot(i).knum,j];
            end
        end
    end
    
    if(knot(i).vertnum==2) % knot sur le bord
        if(knot(i).inum(1)==1 && knot(i).inum(2)==1)
            %vert au sud
            knot(i).inum(1)=0;
            knot(i).inum(2)=0;
        end
        if(knot(i).inum(1)==zone(knot(i).knum(1)).Nx && knot(i).inum(2)==zone(knot(i).knum(2)).Nx)
            %vert au nord
            knot(i).inum(1)=zone(knot(i).knum(1)).Nx+1;
            knot(i).inum(2)=zone(knot(i).knum(2)).Nx+1;
        end
        if(knot(i).jnum(1)==1 && knot(i).jnum(2)==1)
            %vert a l'ouest
            knot(i).jnum(1)=0;
            knot(i).jnum(2)=0;
        end
        if(knot(i).jnum(1)==zone(knot(i).knum(1)).Nz && knot(i).jnum(2)==zone(knot(i).knum(2)).Nz)
            %vert a l'est
            knot(i).jnum(1)=zone(knot(i).knum(1)).Nz+1;
            knot(i).jnum(2)=zone(knot(i).knum(2)).Nz+1;
        end
    end
    
    fprintf(fid,strcat(num2str(i),'\t\t!knot number\n'));
    fprintf(fid,strcat(num2str(knot(i).vertnum),'\t\t!number of vertices\n'));
    for i2=1:knot(i).vertnum
        fprintf(fid,strcat(num2str(knot(i).knum(i2)),'\t'));
    end
    fprintf(fid,'!zone number\n');
    for i2=1:knot(i).vertnum
        fprintf(fid,strcat(num2str(knot(i).inum(i2)),'\t'));
    end
    fprintf(fid,'!line number\n');
    for i2=1:knot(i).vertnum
        fprintf(fid,strcat(num2str(knot(i).jnum(i2)),'\t'));
    end
    fprintf(fid,'!colomn number\n');
end
fclose(fid)

!mkdir triangles
cd triangles
for i=1:16
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    filename=strcat('triA_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(1:Nx,1:Nz);
    fid=fopen(filename,'w');
    for j=1:Nx
        str=sprintf('%d\t',A(j,:));
        fprintf(fid,str);
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    filename=strcat('triB_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(2:Nx+1,1:Nz);
    fid=fopen(filename,'w');
    for j=1:Nx
        str=sprintf('%d\t',A(j,:));
        fprintf(fid,str);
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    filename=strcat('triC_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(2:Nx+1,2:Nz+1);
    fid=fopen(filename,'w');
    for j=1:Nx
        str=sprintf('%d\t',A(j,:));
        fprintf(fid,str);
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    filename=strcat('triD_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(1:Nx,2:Nz+1);
    fid=fopen(filename,'w');
    for j=1:Nx
        str=sprintf('%d\t',A(j,:));
        fprintf(fid,str);
        fprintf(fid,'\n');
    end
    fclose(fid);
end

fid3=fopen('soledge2D.npco_char','w');
for i=1:ntri
    fprintf(fid3,'%d\t',i);
    fprintf(fid3,'%12.7e\t',knot(i).R);
    fprintf(fid3,'%12.7e\n',knot(i).Z);
end
fclose(fid3)

cd ..
