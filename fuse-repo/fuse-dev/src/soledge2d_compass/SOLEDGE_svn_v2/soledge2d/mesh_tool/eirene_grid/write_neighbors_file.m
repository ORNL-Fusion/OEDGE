global eirene
ntriangles_=eirene.ntriangles_;
triangles_=eirene.triangles_;

fid = fopen('soledge2D.neighbors','w');
fprintf(fid,'%d\n',ntriangles_);
for n=1:ntriangles_
    fprintf(fid,'%d\t',n);
    fprintf(fid,'%d\t',triangles_(n).neigh1);
    fprintf(fid,'%d\t',triangles_(n).typeneigh1);
    fprintf(fid,'%d\t',triangles_(n).BC1);
    fprintf(fid,'%d\t',triangles_(n).neigh2);
    fprintf(fid,'%d\t',triangles_(n).typeneigh2);
    fprintf(fid,'%d\t',triangles_(n).BC2);
    fprintf(fid,'%d\t',triangles_(n).neigh3);
    fprintf(fid,'%d\t',triangles_(n).typeneigh3);
    fprintf(fid,'%d\t',triangles_(n).BC3);
    fprintf(fid,'%d\t',0);
    fprintf(fid,'%d\n',0);
end
fclose(fid);