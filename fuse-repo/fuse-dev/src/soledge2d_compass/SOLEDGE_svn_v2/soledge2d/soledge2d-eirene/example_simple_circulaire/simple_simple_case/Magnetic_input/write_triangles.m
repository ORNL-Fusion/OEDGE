fid=fopen('vert2knots.txt','w');
fprintf(fid,strcat(num2str(ntri),'\t\t!number of knots\n'));
for i=1:ntri
    knot(i).vertnum=0;
    knot(i).inum=[];
    knot(i).jnum=[];
    knot(i).knum=[];
    for j=1:2
        Nx=zone(j).Nx;
        Nz=zone(j).Nz;
        
        [m,n]=find(zone(j).ntrinum==i);
        if(not(isempty([m,n])))
            knot(i).R=zone(j).Rtri(m(1),n(1));
            knot(i).Z=zone(j).Ztri(m(1),n(1));
        end
        am=length(m);
        for np=1:am
            %noeud au milieu d'une zone
            if(m(np)>=2)
                if(m(np)<=Nx)
                    if(n(np)>=2)
                        if(n(np)<=Nz)
                            %pt1
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,m(np)-1];
                            knot(i).jnum=[knot(i).jnum,n(np)-1];
                            knot(i).knum=[knot(i).knum,j];
                            %pt2
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,m(np)];
                            knot(i).jnum=[knot(i).jnum,n(np)-1];
                            knot(i).knum=[knot(i).knum,j];
                            %pt3
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,m(np)];
                            knot(i).jnum=[knot(i).jnum,n(np)];
                            knot(i).knum=[knot(i).knum,j];
                            %pt4
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,m(np)-1];
                            knot(i).jnum=[knot(i).jnum,n(np)];
                            knot(i).knum=[knot(i).knum,j];
                        end
                    end
                end
            end
            if(m(np)==1) %knot au sud
                if(n(np)<=Nz)
                    if(n(np)>=2)
                        %pt1
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,1];
                        knot(i).jnum=[knot(i).jnum,n(np)];
                        knot(i).knum=[knot(i).knum,j];
                        %pt2
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,1];
                        knot(i).jnum=[knot(i).jnum,n(np)-1];
                        knot(i).knum=[knot(i).knum,j];
                    end
                end
                if(n(np)==1)
                        %pt1
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,1];
                        knot(i).jnum=[knot(i).jnum,1];
                        knot(i).knum=[knot(i).knum,j];
                end
                if(n(np)==Nz+1)
                    %pt1
                    knot(i).vertnum=knot(i).vertnum+1;
                    knot(i).inum=[knot(i).inum,1];
                    knot(i).jnum=[knot(i).jnum,Nz];
                    knot(i).knum=[knot(i).knum,j];
                end
            end
            
            if(m(np)==Nx+1) %knot au nord
                if(n(np)<=Nz)
                    if(n(np)>=2)
                        %pt1
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,Nx];
                        knot(i).jnum=[knot(i).jnum,n(np)];
                        knot(i).knum=[knot(i).knum,j];
                        %pt2
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,Nx];
                        knot(i).jnum=[knot(i).jnum,n(np)-1];
                        knot(i).knum=[knot(i).knum,j];
                    end
                end
                if(n(np)==1)
                        %pt1
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,Nx];
                        knot(i).jnum=[knot(i).jnum,1];
                        knot(i).knum=[knot(i).knum,j];
                end
                if(n(np)==Nz+1)
                    %pt1
                    knot(i).vertnum=knot(i).vertnum+1;
                    knot(i).inum=[knot(i).inum,Nx];
                    knot(i).jnum=[knot(i).jnum,Nz];
                    knot(i).knum=[knot(i).knum,j];
                end
            end
            
            if(n(np)==1) %knot a l'ouest
                if(m(np)<=Nx)
                    if(m(np)>=2)
                        %pt1
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,m(np)];
                        knot(i).jnum=[knot(i).jnum,1];
                        knot(i).knum=[knot(i).knum,j];
                        %pt2
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,m(np)-1];
                        knot(i).jnum=[knot(i).jnum,1];
                        knot(i).knum=[knot(i).knum,j];
                    end
                end
            end
            
            if(n(np)==Nz+1) %knot a l'est
                if(m(np)<=Nx)
                    if(m(np)>=2)
                        %pt1
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,m(np)];
                        knot(i).jnum=[knot(i).jnum,Nz];
                        knot(i).knum=[knot(i).knum,j];
                        %pt2
                        knot(i).vertnum=knot(i).vertnum+1;
                        knot(i).inum=[knot(i).inum,m(np)-1];
                        knot(i).jnum=[knot(i).jnum,Nz];
                        knot(i).knum=[knot(i).knum,j];
                    end
                end
            end
            
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
for i=1:2
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
fprintf(fid3,'%d\n',ntri);
for i=1:ntri
    fprintf(fid3,'%d\t',i);
    fprintf(fid3,'%12.7e\t',100.*knot(i).R);
    fprintf(fid3,'%12.7e\n',100.*knot(i).Z);
end
fclose(fid3)

cd ..
