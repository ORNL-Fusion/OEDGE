global eirene

nstep=1;
for n=1:ntriangles
    triangles(n).step=0;
end
triangles(floor(ntriangles/2)).step=1;

keep_going=true;
while(keep_going)
    keep_going=false;
    for n=1:ntriangles
        if(triangles(n).step==0)
            if(triangles(n).BC1==0)
                if(triangles(triangles(n).neigh1).step==1)
                    triangles(n).step=1;
                    keep_going=true;
                end
            end
            if(triangles(n).BC2==0)
                if(triangles(triangles(n).neigh2).step==1)
                    triangles(n).step=1;
                    keep_going=true;
                end
            end
            if(triangles(n).BC3==0)
                if(triangles(triangles(n).neigh3).step==1)
                    triangles(n).step=1;
                    keep_going=true;
                end
            end
        end
    end
    for n=ntriangles:-1:1
        if(triangles(n).step==0)
            if(triangles(n).BC1==0)
                if(triangles(triangles(n).neigh1).step==1)
                    triangles(n).step=1;
                    keep_going=true;
                end
            end
            if(triangles(n).BC2==0)
                if(triangles(triangles(n).neigh2).step==1)
                    triangles(n).step=1;
                    keep_going=true;
                end
            end
            if(triangles(n).BC3==0)
                if(triangles(triangles(n).neigh3).step==1)
                    triangles(n).step=1;
                    keep_going=true;
                end
            end
        end
    end


    nstep=nstep+1
end