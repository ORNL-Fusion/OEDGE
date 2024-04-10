fid=fopen('soledge2D.elemente','w');
fprintf(fid,'%d\n',ntriangles_);
for n=1:ntriangles_
    fprintf(fid,'%d\t',n);
    fprintf(fid,'%d\t',triangles_(n).p1);
    fprintf(fid,'%d\t',triangles_(n).p2);
    fprintf(fid,'%d\n',triangles_(n).p3);
end
fclose(fid);
