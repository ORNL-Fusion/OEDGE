close all
figure(1)
hold on
cd ../triangles


fid=fopen('soledge2D.npco_char');
file=fscanf(fid,'%f');
fclose(fid);
ntri_=file(1);
knots=zeros(3,ntri_);
for n=1:ntri_
    knots(n,2)=file((n-1)*3+1+2);
    knots(n,3)=file((n-1)*3+1+3);
end
fclose(fid)

% knots=load('soledge2D.npco_char');

fid=fopen('soledge2D.elemente');
file=fscanf(fid,'%d');
fclose(fid);
ntri_=file(1);
for n=1:ntri_
    tri_(n).k1=file((n-1)*4+1+2);
    tri_(n).k2=file((n-1)*4+1+3);
    tri_(n).k3=file((n-1)*4+1+4);
end


% fid=fopen('soledge2D.neighbors');
% file=fscanf(fid,'%d');
% fclose(fid);
Neigh=load('soledge2D.neighbors');
triV1=Neigh(:,4);
triV2=Neigh(:,7);
triV3=Neigh(:,10);

for n=1:ntri_
%     axis equal
    if(triV1(n)==0)
        plot([knots(tri_(n).k1,2),knots(tri_(n).k2,2)],[knots(tri_(n).k1,3),knots(tri_(n).k2,3)],'k-')
    elseif(triV1(n)==1)
        plot([knots(tri_(n).k1,2),knots(tri_(n).k2,2)],[knots(tri_(n).k1,3),knots(tri_(n).k2,3)],'r-','LineWidth',5)
    elseif(triV1(n)==2)
        plot([knots(tri_(n).k1,2),knots(tri_(n).k2,2)],[knots(tri_(n).k1,3),knots(tri_(n).k2,3)],'g-','LineWidth',5)
    end
    
    if(triV2(n)==0)
        plot([knots(tri_(n).k2,2),knots(tri_(n).k3,2)],[knots(tri_(n).k2,3),knots(tri_(n).k3,3)],'k-')
    elseif(triV2(n)==1)
        plot([knots(tri_(n).k2,2),knots(tri_(n).k3,2)],[knots(tri_(n).k2,3),knots(tri_(n).k3,3)],'r-','LineWidth',5)
    elseif(triV2(n)==2)
        plot([knots(tri_(n).k2,2),knots(tri_(n).k3,2)],[knots(tri_(n).k2,3),knots(tri_(n).k3,3)],'g-','LineWidth',5)
    end
    
    if(triV3(n)==0)
        plot([knots(tri_(n).k3,2),knots(tri_(n).k1,2)],[knots(tri_(n).k3,3),knots(tri_(n).k1,3)],'k-')
    elseif(triV3(n)==1)
        plot([knots(tri_(n).k3,2),knots(tri_(n).k1,2)],[knots(tri_(n).k3,3),knots(tri_(n).k1,3)],'r-','LineWidth',5)
    elseif(triV3(n)==2)
        plot([knots(tri_(n).k3,2),knots(tri_(n).k1,2)],[knots(tri_(n).k3,3),knots(tri_(n).k1,3)],'g-','LineWidth',5)
    end

    if(mod(n,100)==0)
    drawnow
    end
    
end

cd ../Magnetic_input
