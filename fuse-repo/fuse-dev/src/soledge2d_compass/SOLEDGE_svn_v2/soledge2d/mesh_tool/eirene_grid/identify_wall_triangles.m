global eirene;

wall_tri=[];
for n=1:ntriangles_
   if((triangles_(n).BC1==3)||(triangles_(n).BC2==3)||(triangles_(n).BC3==3))
      wall_tri=[wall_tri;n]; 
   end
end


ntriwall=length(wall_tri);
fid=fopen('soledge2D.backinterp','w');
fprintf(fid,'%d\n',ntriwall);
eirene.wall.num=ntriwall;
eirene.wall.tri=wall_tri;
eirene.wall.tri_i=zeros(size(wall_tri));
eirene.wall.tri_j=zeros(size(wall_tri));
eirene.wall.tri_k=zeros(size(wall_tri));

for n=1:ntriwall
    eirene.wall.tri_k(n)=triangles_(wall_tri(n)).k;
    eirene.wall.tri_i(n)=triangles_(wall_tri(n)).i;
    eirene.wall.tri_j(n)=triangles_(wall_tri(n)).j;
    fprintf(fid,'%d\t',wall_tri(n));
    fprintf(fid,'%d\t',triangles_(wall_tri(n)).k);
    fprintf(fid,'%d\t',triangles_(wall_tri(n)).i);
    fprintf(fid,'%d\n',triangles_(wall_tri(n)).j);
end
fclose(fid)

fid=fopen('soledge2D.backinterp_full','w');
for n=1:ntriangles_
    fprintf(fid,'%d\t',n);
    fprintf(fid,'%d\t',triangles_(n).k);
    fprintf(fid,'%d\t',triangles_(n).i);
    fprintf(fid,'%d\n',triangles_(n).j);
end
fclose(fid);

