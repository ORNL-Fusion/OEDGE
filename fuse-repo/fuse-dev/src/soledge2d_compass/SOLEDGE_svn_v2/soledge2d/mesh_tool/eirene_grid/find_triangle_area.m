global megazone

L=[];
for k=1:megazone.num
    if(megazone.mz(k).isperiodic)
        L=[L,megazone.mz(k).list];
    end
end

for n=1:ntriangles
    if(sum(find(triangles(n).k==L))~=0)
        triangles(n).area=1; %in the core
    else
        triangles(n).area=0; %at the periphery
    end
end
