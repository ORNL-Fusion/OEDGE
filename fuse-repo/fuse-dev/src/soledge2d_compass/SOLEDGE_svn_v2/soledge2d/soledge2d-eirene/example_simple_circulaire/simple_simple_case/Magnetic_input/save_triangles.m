ntriangle=0;

tri_knots=[];

%attributing triangle numbers and vertices
for i=1:2
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    zone(i).triangles=zeros(Nx,Nz,2);
    zone(i).trivert=zeros(Nx,Nz,2,3);
    zone(i).outfluxtri=zeros(Nx,Nz,4);
    zone(i).outfluxtrivert=zeros(Nx,Nz,4);
    for k=1:Nx
        for j=1:Nz
            %premier triangle
            nt1=zone(i).ntrinum(k,j);
            nt2=zone(i).ntrinum(k+1,j);
            nt3=zone(i).ntrinum(k+1,j+1);
            ntriangle=ntriangle+1;
            tri_knots=[tri_knots;nt1,nt2,nt3];
            
            zone(i).triangles(k,j,1)=ntriangle;
            zone(i).trivert(k,j,1,1)=nt1;
            zone(i).trivert(k,j,1,2)=nt2;
            zone(i).trivert(k,j,1,3)=nt3;
            
            zone(i).outfluxtri(k,j,1)=ntriangle;
            zone(i).outfluxtrivert(k,j,1)=2;
            zone(i).outfluxtri(k,j,4)=ntriangle;
            zone(i).outfluxtrivert(k,j,4)=3;
             
            %deuxieme triangle
            nt1=zone(i).ntrinum(k,j);
            nt2=zone(i).ntrinum(k+1,j+1);
            nt3=zone(i).ntrinum(k,j+1);
            ntriangle=ntriangle+1
            tri_knots=[tri_knots;nt1,nt2,nt3];
            
            zone(i).triangles(k,j,2)=ntriangle;
            zone(i).trivert(k,j,2,1)=nt1;
            zone(i).trivert(k,j,2,2)=nt2;
            zone(i).trivert(k,j,2,3)=nt3;
            
            zone(i).outfluxtri(k,j,2)=ntriangle;
            zone(i).outfluxtrivert(k,j,2)=1;
            zone(i).outfluxtri(k,j,3)=ntriangle;
            zone(i).outfluxtrivert(k,j,3)=2;
        end
    end
end

%attributing neighbors
cd ..
fid=fopen('Neighbors.txt');
fid2=fopen('MagNeighbors.txt');
cd Magnetic_input
for i=1:2
   Nx=zone(i).Nx;
   Nz=zone(i).Nz;
   A=fscanf(fid,'%d',4);
   B=fscanf(fid2,'%d',4);
   % neighboring triangle on side 1,2,3 of triangle 1,2
   zone(i).neightri=zeros(Nx,Nz,2,3);
   % entrance side on neighboring triangle on side 1,2,3 for triangle 1,2
   zone(i).neighsid=zeros(Nx,Nz,2,3);
   % surface model on side 1,2,3 of triangle 1,2; 0 means no material
   % surface
   zone(i).neightyp=zeros(Nx,Nz,2,3);
    
   for k=1:Nx
       for j=1:Nz
           % k-1,j neighbor
           if(k==1)
                if(B(2)==1) %BC au sud 
                   zone(i).neightri(k,j,2,3)=0;
                   zone(i).neighsid(k,j,2,3)=0;
                   zone(i).neightyp(k,j,2,3)=1;
                else
                   zone(i).neightri(k,j,2,3)=zone(A(2)).triangles(zone(A(2)).Nx,j,1);
                   zone(i).neighsid(k,j,2,3)=2;
                   zone(i).neightyp(k,j,2,3)=0;
                end
            else
                zone(i).neightri(k,j,2,3)=zone(i).triangles(k-1,j,1);
                zone(i).neighsid(k,j,2,3)=2;
                zone(i).neightyp(k,j,2,3)=0;
           end
            
            % k,j-1 neighbor
            if(j==1) 
               if(B(4)==1) %BC à l'ouest
                   zone(i).neightri(k,j,1,1)=0;
                   zone(i).neighsid(k,j,1,1)=0;
                   zone(i).neightyp(k,j,1,1)=1;
               else
                   zone(i).neightri(k,j,1,1)=zone(A(4)).triangles(k,zone(A(4)).Nz,2);
                   zone(i).neighsid(k,j,1,1)=2;
                   zone(i).neightyp(k,j,1,1)=0;
               end
            else
                zone(i).neightri(k,j,1,1)=zone(i).triangles(k,j-1,2);
                zone(i).neighsid(k,j,1,1)=2;
                zone(i).neightyp(k,j,1,1)=0;
            end
            
           % k+1,j neighbor
            if(k==Nx)
               if(B(1)==1) %BC au nord
                  zone(i).neightri(k,j,1,2)=0;
                  zone(i).neighsid(k,j,1,2)=0;
                  zone(i).neightyp(k,j,1,2)=1;
               else
                  zone(i).neightri(k,j,1,2)=zone(A(1)).triangles(1,j,2);
                  zone(i).neighsid(k,j,1,2)=3;
                  zone(i).neightyp(k,j,1,2)=0; 
               end
            else
                zone(i).neightri(k,j,1,2)=zone(i).triangles(k+1,j,2);
                zone(i).neighsid(k,j,1,2)=3;
                zone(i).neightyp(k,j,1,2)=0;
            end
            % k, j+1 neighbor
            if(j==Nz)
                if(B(3)==1) 
                    zone(i).neightri(k,j,2,2)=0;
                    zone(i).neighsid(k,j,2,2)=0;
                    zone(i).neightyp(k,j,2,2)=1;
                else
                    zone(i).neightri(k,j,2,2)=zone(A(3)).triangles(k,1,1);
                    zone(i).neighsid(k,j,2,2)=1;
                    zone(i).neightyp(k,j,2,2)=0;
                end
            else
                zone(i).neightri(k,j,2,2)=zone(i).triangles(k,j+1,1);
                zone(i).neighsid(k,j,2,2)=1;
                zone(i).neightyp(k,j,2,2)=0;
            end
            % this is the boundary inside a quadrangle
            zone(i).neightri(k,j,1,3)=zone(i).triangles(k,j,2);
            zone(i).neighsid(k,j,1,3)=1;
            zone(i).neightri(k,j,2,1)=zone(i).triangles(k,j,1);
            zone(i).neighsid(k,j,2,1)=3;
       end
   end
end
   fclose(fid);
   fclose(fid2);
   
   %saving
%    SV=[];
cd triangles
   fid=fopen('soledge2D.neighbors','w');
   fid2=fopen('soledge2D.elemente','w');
   fid4=fopen('soledge2D.zones','w');
   fprintf(fid4,'%d\n',ntriangle);
   fprintf(fid,'%s\n',' ');
   fprintf(fid2,'%d\n',ntriangle);
   
   for i=1:2
    Nx=zone(i).Nx;
    Nz=zone(i).Nz;
    for k=1:Nx
        for j=1:Nz
            %on ajoute le premier triangle
%             SV=[SV;[zone(i).triangles(k,j,1),zone(i).neightri(k,j,1,1),3,zone(i).neightype(k,j,1,1)...
%                 zone(i).neightri(k,j,1,2),2,zone(i).neightype(k,j,1,2),...
%                 zone(i).neightri(k,j,1,3),1,zone(i).neightype(k,j,1,3),0,0]];
            
            fprintf(fid4,'%d\t',zone(i).triangles(k,j,1));
            fprintf(fid4,'%d\t',i);
            fprintf(fid4,'%d\t',k);
            fprintf(fid4,'%d\t',j);
            fprintf(fid4,'%d\t',4);
            fprintf(fid4,'%d\t',1);
            fprintf(fid4,'%d\n',0);
            
            fprintf(fid,'%d\t',zone(i).triangles(k,j,1));
            %premier voisin
            fprintf(fid,'%d\t',zone(i).neightri(k,j,1,1));
            fprintf(fid,'%d\t',zone(i).neighsid(k,j,1,1));
            fprintf(fid,'%d\t',zone(i).neightyp(k,j,1,1));
            %deuxieme voisin
            fprintf(fid,'%d\t',zone(i).neightri(k,j,1,2));
            fprintf(fid,'%d\t',zone(i).neighsid(k,j,1,2));
            fprintf(fid,'%d\t',zone(i).neightyp(k,j,1,2));
            %troisieme voisin
            fprintf(fid,'%d\t',zone(i).neightri(k,j,1,3));
            fprintf(fid,'%d\t',zone(i).neighsid(k,j,1,3));
            fprintf(fid,'%d\t',zone(i).neightyp(k,j,1,3));
            %coordonnée ???
            fprintf(fid,'%d\t',0);
            fprintf(fid,'%d\n',0);
            
            fprintf(fid2,'%d\t',zone(i).triangles(k,j,1));
            fprintf(fid2,'%d\t',zone(i).trivert(k,j,1,1));
            fprintf(fid2,'%d\t',zone(i).trivert(k,j,1,2));
            fprintf(fid2,'%d\n',zone(i).trivert(k,j,1,3));
            
            %check trigo
            vector1x=knot(zone(i).trivert(k,j,1,2)).R-knot(zone(i).trivert(k,j,1,1)).R;
            vector1y=knot(zone(i).trivert(k,j,1,2)).Z-knot(zone(i).trivert(k,j,1,1)).Z;
            
            vector2x=knot(zone(i).trivert(k,j,1,3)).R-knot(zone(i).trivert(k,j,1,1)).R;
            vector2y=knot(zone(i).trivert(k,j,1,3)).Z-knot(zone(i).trivert(k,j,1,1)).Z;
            
            vecprod=vector1x*vector2y-vector1y*vector2x;
            if(vecprod<=0)
               disp(['wrong orientation, triangle',num2str(zone(i).triangles(k,j,1))]) ;
            end
            
            %on ajoute le deuxieme triangle
%             SV=[SV;[zone(i).triangles(k,j,1),zone(i).neightri(k,j,1,1),3,zone(i).neightype(k,j,1,1)...
%                 zone(i).neightri(k,j,1,2),2,zone(i).neightype(k,j,1,2),...
%                 zone(i).neightri(k,j,1,3),1,zone(i).neightype(k,j,1,3),0,0]];
            fprintf(fid,'%d\t',zone(i).triangles(k,j,2));
            fprintf(fid4,'%d\t',zone(i).triangles(k,j,2));
            fprintf(fid4,'%d\t',i);
            fprintf(fid4,'%d\t',k);
            fprintf(fid4,'%d\t',j);
            fprintf(fid4,'%d\t',0);
            fprintf(fid4,'%d\t',3);
            fprintf(fid4,'%d\n',2);
            %premier voisin
            fprintf(fid,'%d\t',zone(i).neightri(k,j,2,1));
            fprintf(fid,'%d\t',zone(i).neighsid(k,j,2,1));
            fprintf(fid,'%d\t',zone(i).neightyp(k,j,2,1));
            %deuxieme voisin
            fprintf(fid,'%d\t',zone(i).neightri(k,j,2,2));
            fprintf(fid,'%d\t',zone(i).neighsid(k,j,2,2));
            fprintf(fid,'%d\t',zone(i).neightyp(k,j,2,2));
            %troisieme voisin
            fprintf(fid,'%d\t',zone(i).neightri(k,j,2,3));
            fprintf(fid,'%d\t',zone(i).neighsid(k,j,2,3));
            fprintf(fid,'%d\t',zone(i).neightyp(k,j,2,3));
            %coordonnée ???
            fprintf(fid,'%d\t',0);
            fprintf(fid,'%d\n',0);
            
            fprintf(fid2,'%d\t',zone(i).triangles(k,j,2));
            fprintf(fid2,'%d\t',zone(i).trivert(k,j,2,1));
            fprintf(fid2,'%d\t',zone(i).trivert(k,j,2,2));
            fprintf(fid2,'%d\n',zone(i).trivert(k,j,2,3));
            
            %check trigo
            vector1x=knot(zone(i).trivert(k,j,2,2)).R-knot(zone(i).trivert(k,j,2,1)).R;
            vector1y=knot(zone(i).trivert(k,j,2,2)).Z-knot(zone(i).trivert(k,j,2,1)).Z;
            
            vector2x=knot(zone(i).trivert(k,j,2,3)).R-knot(zone(i).trivert(k,j,2,1)).R;
            vector2y=knot(zone(i).trivert(k,j,2,3)).Z-knot(zone(i).trivert(k,j,2,1)).Z;
            
            vecprod=vector1x*vector2y-vector1y*vector2x;
            if(vecprod<=0)
               disp(['wrong orientation, triangle',num2str(zone(i).triangles(k,j,2))]) ;
            end
            
        end
    end
   end
fclose(fid)
fclose(fid2)
fclose(fid4)

% for i=1:16
%    %au nord: triangle 2, face 1
%    filename=strcat('out_triN_',num2str(i,'%.3d'));
%    filename=strcat(filename,'.txt');
%    A=zone(i).outfluxtri(:,:,2);
%    fid=fopen(filename,'w');
%    for k=1:zone(i).Nx
%       for j=1:zone(i).Nz
%          fprintf(fid,'%d\t',A(k,j));
%       end
%       fprintf(fid,'\n');
%    end
%    fclose(fid);
%    %au sud: triangle 1, face 3
%    filename=strcat('out_triS_',num2str(i,'%.3d'));
%    filename=strcat(filename,'.txt');
%    A=zone(i).outfluxtri(:,:,4);
%    fid=fopen(filename,'w');
%    for k=1:zone(i).Nx
%       for j=1:zone(i).Nz
%          fprintf(fid,'%d\t',A(k,j));
%       end
%       fprintf(fid,'\n');
%    end
%    fclose(fid);
%    %a l'ouest: triangle 1, face 1
%    filename=strcat('out_triW_',num2str(i,'%.3d'));
%    filename=strcat(filename,'.txt');
%    A=zone(i).outfluxtri(:,:,1);
%    fid=fopen(filename,'w');
%    for k=1:zone(i).Nx
%       for j=1:zone(i).Nz
%          fprintf(fid,'%d\t',A(k,j));
%       end
%       fprintf(fid,'\n');
%    end
%    fclose(fid);
%    %a l'est: triangle 2, face 3
%    filename=strcat('out_triE_',num2str(i,'%.3d'));
%    filename=strcat(filename,'.txt');
%    A=zone(i).outfluxtri(:,:,3);
%    fid=fopen(filename,'w');
%    for k=1:zone(i).Nx
%       for j=1:zone(i).Nz
%          fprintf(fid,'%d\t',A(k,j));
%       end
%       fprintf(fid,'\n');
%    end
%    fclose(fid);
% end

fid=fopen('soledge2D.knots_info','w');
for nt=1:ntri
    [m,p]=find(tri_knots==nt);
    fprintf(fid,'%d\t',nt);
    fprintf(fid,'%d\n',length(m));
    for i=1:length(m)
       fprintf(fid,'%d\t',m(i)); 
    end
    fprintf(fid,'\n');
end
fclose(fid);
cd ..