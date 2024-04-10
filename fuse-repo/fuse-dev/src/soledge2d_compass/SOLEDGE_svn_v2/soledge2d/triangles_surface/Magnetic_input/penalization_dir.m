ntri_new=0;
nknot_new=0;
allowed_knot_list=[];

cd ../Mesh/
for k=1:nzones
    name=strcat('chi_',num2str(k,'%.3d'));
    name=strcat(name,'.txt');
    zone(k).chi=load(name);
end
cd ../Magnetic_input/


for n=1:ntriangle
    if(zone(tri_info(n).k).chi(tri_info(n).i,tri_info(n).j)==0)
        ntri_new=ntri_new+1;
        tri_info(n).newn=ntri_new;
        allowed_knot_list=[allowed_knot_list,tri_info(n).p1,tri_info(n).p2,tri_info(n).p3];
    else
        tri_info(n).newn=0;
    end
end

for n=1:ntri
    if(length(find(allowed_knot_list==n))>0)
        nknot_new=nknot_new+1;
        knot_info(n).newn=nknot_new;
    else
        knot_info(n).newn=0;
    end
end


% copy of write_triangles


cd ../triangles
fid=fopen('vert2knots.txt','w');
fprintf(fid,strcat(num2str(nknot_new),'\t\t!number of knots\n'));
for i=1:ntri
    if(knot_info(i).newn>0)
        knot(i).vertnum=0;
        knot(i).inum=[];
        knot(i).jnum=[];
        knot(i).knum=[];
        for j=1:nzones
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
                                if(zone(j).chi(m(np)-1,n(np)-1)==0)
                                    knot(i).vertnum=knot(i).vertnum+1;
                                    knot(i).inum=[knot(i).inum,m(np)-1];
                                    knot(i).jnum=[knot(i).jnum,n(np)-1];
                                    knot(i).knum=[knot(i).knum,j];
                                end
                                %pt2
                                if(zone(j).chi(m(np),n(np)-1)==0)
                                    knot(i).vertnum=knot(i).vertnum+1;
                                    knot(i).inum=[knot(i).inum,m(np)];
                                    knot(i).jnum=[knot(i).jnum,n(np)-1];
                                    knot(i).knum=[knot(i).knum,j];
                                end
                                %pt3
                                if(zone(j).chi(m(np),n(np))==0)
                                    knot(i).vertnum=knot(i).vertnum+1;
                                    knot(i).inum=[knot(i).inum,m(np)];
                                    knot(i).jnum=[knot(i).jnum,n(np)];
                                    knot(i).knum=[knot(i).knum,j];
                                end
                                %pt4
                                if(zone(j).chi(m(np)-1,n(np))==0)
                                    knot(i).vertnum=knot(i).vertnum+1;
                                    knot(i).inum=[knot(i).inum,m(np)-1];
                                    knot(i).jnum=[knot(i).jnum,n(np)];
                                    knot(i).knum=[knot(i).knum,j];
                                end
                            end
                        end
                    end
                end
                if(m(np)==1) %knot au sud
                    if(n(np)<=Nz)
                        if(n(np)>=2)
                            %pt1
                            if(zone(j).chi(1,n(np))==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,1];
                                knot(i).jnum=[knot(i).jnum,n(np)];
                                knot(i).knum=[knot(i).knum,j];
                            end
                            %pt2
                            if(zone(j).chi(1,n(np)-1)==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,1];
                                knot(i).jnum=[knot(i).jnum,n(np)-1];
                                knot(i).knum=[knot(i).knum,j];
                            end
                        end
                    end
                    if(n(np)==1)
                        %pt1
                        if(zone(j).chi(1,1)==0)
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,1];
                            knot(i).jnum=[knot(i).jnum,1];
                            knot(i).knum=[knot(i).knum,j];
                        end
                    end
                    if(n(np)==Nz+1)
                        %pt1
                        if(zone(j).chi(1,Nz)==0)
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,1];
                            knot(i).jnum=[knot(i).jnum,Nz];
                            knot(i).knum=[knot(i).knum,j];
                        end
                    end
                end
                
                if(m(np)==Nx+1) %knot au nord
                    if(n(np)<=Nz)
                        if(n(np)>=2)
                            %pt1
                            if(zone(j).chi(Nx,n(np))==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,Nx];
                                knot(i).jnum=[knot(i).jnum,n(np)];
                                knot(i).knum=[knot(i).knum,j];
                            end
                            %pt2
                            if(zone(j).chi(Nx,n(np)-1)==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,Nx];
                                knot(i).jnum=[knot(i).jnum,n(np)-1];
                                knot(i).knum=[knot(i).knum,j];
                            end
                        end
                    end
                    if(n(np)==1)
                        %pt1
                        if(zone(j).chi(Nx,1)==0)
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,Nx];
                            knot(i).jnum=[knot(i).jnum,1];
                            knot(i).knum=[knot(i).knum,j];
                        end
                    end
                    if(n(np)==Nz+1)
                        %pt1
                        if(zone(j).chi(Nx,Nz)==0)
                            knot(i).vertnum=knot(i).vertnum+1;
                            knot(i).inum=[knot(i).inum,Nx];
                            knot(i).jnum=[knot(i).jnum,Nz];
                            knot(i).knum=[knot(i).knum,j];
                        end
                    end
                end
                
                if(n(np)==1) %knot a l'ouest
                    if(m(np)<=Nx)
                        if(m(np)>=2)
                            %pt1
                            if(zone(j).chi(m(np),1)==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,m(np)];
                                knot(i).jnum=[knot(i).jnum,1];
                                knot(i).knum=[knot(i).knum,j];
                            end
                            %pt2
                            if(zone(j).chi(m(np)-1,1)==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,m(np)-1];
                                knot(i).jnum=[knot(i).jnum,1];
                                knot(i).knum=[knot(i).knum,j];
                            end
                        end
                    end
                end
                
                if(n(np)==Nz+1) %knot a l'est
                    if(m(np)<=Nx)
                        if(m(np)>=2)
                            %pt1
                            if(zone(j).chi(m(np),Nz)==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,m(np)];
                                knot(i).jnum=[knot(i).jnum,Nz];
                                knot(i).knum=[knot(i).knum,j];
                            end
                            %pt2
                            if(zone(j).chi(m(np)-1,Nz)==0)
                                knot(i).vertnum=knot(i).vertnum+1;
                                knot(i).inum=[knot(i).inum,m(np)-1];
                                knot(i).jnum=[knot(i).jnum,Nz];
                                knot(i).knum=[knot(i).knum,j];
                            end
                        end
                    end
                end
                
            end
        end
        
        fprintf(fid,strcat(num2str(knot_info(i).newn),'\t\t!knot number\n'));
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
end
fclose(fid)

for i=1:nzones
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    filename=strcat('triA_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(1:Nx,1:Nz);
    fid=fopen(filename,'w');
    for j=1:Nx
        for k=1:Nz
            if(zone(i).chi(j,k)==0)
                fprintf(fid,'%d\t',knot_info(A(j,k)).newn);
            else
                fprintf(fid,'%d\t',0);
            end
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    filename=strcat('triB_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(2:Nx+1,1:Nz);
    fid=fopen(filename,'w');
    for j=1:Nx
        for k=1:Nz
            if(zone(i).chi(j,k)==0)
                fprintf(fid,'%d\t',knot_info(A(j,k)).newn);
            else
                fprintf(fid,'%d\t',0);
            end
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    filename=strcat('triC_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(2:Nx+1,2:Nz+1);
    fid=fopen(filename,'w');
    for j=1:Nx
        for k=1:Nz
            if(zone(i).chi(j,k)==0)
                fprintf(fid,'%d\t',knot_info(A(j,k)).newn);
            else
                fprintf(fid,'%d\t',0);
            end
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    filename=strcat('triD_',num2str(i,'%.3d'));
    filename=strcat(filename,'.txt');
    A=zone(i).ntrinum(1:Nx,2:Nz+1);
    fid=fopen(filename,'w');
    for j=1:Nx
        for k=1:Nz
            if(zone(i).chi(j,k)==0)
                fprintf(fid,'%d\t',knot_info(A(j,k)).newn);
            else
                fprintf(fid,'%d\t',0);
            end
        end
        fprintf(fid,'\n');
    end
    fclose(fid);
end

fid3=fopen('soledge2D.npco_char','w');
fprintf(fid3,'%d\n',nknot_new);
for i=1:ntri
    if(knot_info(i).newn>0)
        fprintf(fid3,'%d\t',knot_info(i).newn);
        fprintf(fid3,'%12.7e\t',knot(i).R.*100);
        fprintf(fid3,'%12.7e\n',knot(i).Z.*100);
    end
end
fclose(fid3)

% copy of save_triangles

fid=fopen('soledge2D.neighbors','w');
fid2=fopen('soledge2D.elemente','w');
fid4=fopen('soledge2D.zones','w');
fprintf(fid4,'%d\n',ntri_new);
fprintf(fid,'%s\n',' ');
fprintf(fid2,'%d\n',ntri_new);


for i=1:nzones
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    for k=1:Nx
        for j=1:Nz
            %on ajoute le premier triangle
            if(tri_info(zone(i).triangles(k,j,1)).newn>0)
                fprintf(fid4,'%d\t',tri_info(zone(i).triangles(k,j,1)).newn);
                fprintf(fid4,'%d\t',i);
                fprintf(fid4,'%d\t',k);
                fprintf(fid4,'%d\t',j);
                fprintf(fid4,'%d\t',4);
                fprintf(fid4,'%d\t',1);
                fprintf(fid4,'%d\n',0);
                
                fprintf(fid,'%d\t',tri_info(zone(i).triangles(k,j,1)).newn);
                %premier voisin
                if(zone(i).neightri(k,j,1,1)>0)
                    if(tri_info(zone(i).neightri(k,j,1,1)).newn>0)
                        fprintf(fid,'%d\t',tri_info(zone(i).neightri(k,j,1,1)).newn);
                        fprintf(fid,'%d\t',zone(i).neighsid(k,j,1,1));
                        fprintf(fid,'%d\t',zone(i).neightyp(k,j,1,1));
                    else
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',1);
                    end
                else
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',2);
                end
                %deuxieme voisin
                if(zone(i).neightri(k,j,1,2)>0)
                    if(tri_info(zone(i).neightri(k,j,1,2)).newn>0)
                        fprintf(fid,'%d\t',tri_info(zone(i).neightri(k,j,1,2)).newn);
                        fprintf(fid,'%d\t',zone(i).neighsid(k,j,1,2));
                        fprintf(fid,'%d\t',zone(i).neightyp(k,j,1,2));
                    else
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',1);
                    end
                else
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',2);
                end
                %troisieme voisin
                if(zone(i).neightri(k,j,1,3)>0)
                    if(tri_info(zone(i).neightri(k,j,1,3)).newn>0)
                        fprintf(fid,'%d\t',tri_info(zone(i).neightri(k,j,1,3)).newn);
                        fprintf(fid,'%d\t',zone(i).neighsid(k,j,1,3));
                        fprintf(fid,'%d\t',zone(i).neightyp(k,j,1,3));
                    else
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',1);
                    end
                else
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',2);
                end
                %coordonnée ???
                fprintf(fid,'%d\t',0);
                fprintf(fid,'%d\n',0);
                
                fprintf(fid2,'%d\t',tri_info(zone(i).triangles(k,j,1)).newn);
                fprintf(fid2,'%d\t',knot_info(zone(i).trivert(k,j,1,1)).newn);
                fprintf(fid2,'%d\t',knot_info(zone(i).trivert(k,j,1,2)).newn);
                fprintf(fid2,'%d\n',knot_info(zone(i).trivert(k,j,1,3)).newn);
                
            end
            
            
            %on ajoute le deuxieme triangle
            if(tri_info(zone(i).triangles(k,j,2)).newn>0)
                fprintf(fid4,'%d\t',tri_info(zone(i).triangles(k,j,2)).newn);
                fprintf(fid4,'%d\t',i);
                fprintf(fid4,'%d\t',k);
                fprintf(fid4,'%d\t',j);
                fprintf(fid4,'%d\t',0);
                fprintf(fid4,'%d\t',3);
                fprintf(fid4,'%d\n',2);
                
                fprintf(fid,'%d\t',tri_info(zone(i).triangles(k,j,2)).newn);
                %premier voisin
                if(zone(i).neightri(k,j,2,1)>0)
                    if(tri_info(zone(i).neightri(k,j,2,1)).newn>0)
                        fprintf(fid,'%d\t',tri_info(zone(i).neightri(k,j,2,1)).newn);
                        fprintf(fid,'%d\t',zone(i).neighsid(k,j,2,1));
                        fprintf(fid,'%d\t',zone(i).neightyp(k,j,2,1));
                    else
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',1);
                    end
                else
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',2);
                end
                %deuxieme voisin
                if(zone(i).neightri(k,j,2,2)>0)
                    if(tri_info(zone(i).neightri(k,j,2,2)).newn>0)
                        fprintf(fid,'%d\t',tri_info(zone(i).neightri(k,j,2,2)).newn);
                        fprintf(fid,'%d\t',zone(i).neighsid(k,j,2,2));
                        fprintf(fid,'%d\t',zone(i).neightyp(k,j,2,2));
                    else
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',1);
                    end
                else
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',2);
                end
                %troisieme voisin
                if(zone(i).neightri(k,j,2,3)>0)
                    if(tri_info(zone(i).neightri(k,j,2,3)).newn>0)
                        fprintf(fid,'%d\t',tri_info(zone(i).neightri(k,j,2,3)).newn);
                        fprintf(fid,'%d\t',zone(i).neighsid(k,j,2,3));
                        fprintf(fid,'%d\t',zone(i).neightyp(k,j,2,3));
                    else
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',0);
                        fprintf(fid,'%d\t',1);
                    end
                else
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',0);
                    fprintf(fid,'%d\t',2);
                end
                %coordonnée ???
                fprintf(fid,'%d\t',0);
                fprintf(fid,'%d\n',0);
                
                fprintf(fid2,'%d\t',tri_info(zone(i).triangles(k,j,2)).newn);
                fprintf(fid2,'%d\t',knot_info(zone(i).trivert(k,j,2,1)).newn);
                fprintf(fid2,'%d\t',knot_info(zone(i).trivert(k,j,2,2)).newn);
                fprintf(fid2,'%d\n',knot_info(zone(i).trivert(k,j,2,3)).newn);
            end
        end
    end
end
fclose(fid)
fclose(fid2)
fclose(fid4)

tri_knots_new=zeros(ntri_new,3);
for n=1:ntriangle
    if(tri_info(n).newn>0)
        tri_knots_new(tri_info(n).newn,1)=knot_info(tri_knots(n,1)).newn;
        tri_knots_new(tri_info(n).newn,2)=knot_info(tri_knots(n,2)).newn;
        tri_knots_new(tri_info(n).newn,3)=knot_info(tri_knots(n,3)).newn;
    end
end


fid=fopen('soledge2D.knots_info','w');
for nt=1:nknot_new
    [m,p]=find(tri_knots_new==nt);
    fprintf(fid,'%d\t',nt);
    fprintf(fid,'%d\n',length(m));
    for i=1:length(m)
       fprintf(fid,'%d\t',m(i)); 
    end
    fprintf(fid,'\n');
end
fclose(fid);

cd ../Magnetic_input
