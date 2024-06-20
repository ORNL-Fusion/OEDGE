fid=fopen('soledge2D.npco_char','w');
fprintf(fid,'%d\n',nknots_);
for n=1:nknots_
   fprintf(fid,'%d\t',n);
   fprintf(fid,'%12.7e\t',100*R_(n));
   fprintf(fid,'%12.7e\n',100*Z_(n));
end
fclose(fid);
