clear zones
global zones

zones.num=0;


for k=1:ncotes
    
    p1=cotes(k).p1.coord;
    p2=cotes(k).p2.coord;
    list1=[];
    list2=[];
    for k2=1:ncotes
        if(k2~=k)
            comp1=[p1(2),p1(3),p1(4),p2(2),p2(3),p2(4)];
            p1_=cotes(k2).p1.coord;
            p2_=cotes(k2).p2.coord;
            comp2=[p1_(2),p1_(3),p1_(4),p2_(2),p2_(3),p2_(4)];
            somme=0;
%             if(comp1(2)==comp1(5)) %periodic segment detected
%                 for k3=1:6
%                     if((k3~=2)&&(k3~=5))
%                         if(comp2(k3)==comp1(k3))
%                             somme=somme+1;
%                         end
%                     end
%                 end
%                 if(comp2(2)==comp2(5))
%                     if(comp2(2)==comp1(2))
%                         somme=somme+2;
%                     end
%                 end
%                 if((comp2(5)==-1)||(comp2(2)==-1))
%                     somme=somme+2;
%                 end
%             else
%                 if(comp2(2)==comp2(5)) %periodic segment detected
%                     for k3=1:6
%                         if((k3~=2)&&(k3~=5))
%                             if(comp1(k3)==comp2(k3))
%                                 somme=somme+1;
%                             end
%                         end
%                     end
%                     if(comp1(2)==comp1(5))
%                         if(comp1(2)==comp2(2))
%                             somme=somme+2;
%                         end
%                     end
%                     if((comp1(5)==-1)||(comp1(2)==-1))
%                         somme=somme+2;
%                     end
%                 else
                    for k3=1:6
                        if(comp1(k3)==comp2(k3))
                            somme=somme+1;
                        else
                            if((k3==2)||(k3==5))
                                if((comp1(k3)==-1)||(comp2(k3)==-1))
                                    somme=somme+1;
                                    %joker
                                end
                            end
                        end
                    end
%                 end
%             end
            if(somme==6)
                list1=[list1,k2];
            end
            comp2=[p2_(2),p2_(3),p2_(4),p1_(2),p1_(3),p1_(4)];
            somme=0;
%             if(comp1(2)==comp1(5)) %periodic segment detected
%                 for k3=1:6
%                     if((k3~=2)&&(k3~=5))
%                         if(comp2(k3)==comp1(k3))
%                             somme=somme+1;
%                         end
%                     end
%                 end
%                 if(comp2(2)==comp2(5))
%                     if(comp2(2)==comp1(2))
%                         somme=somme+2;
%                     end
%                 end
%                 if((comp2(5)==-1)||(comp2(2)==-1))
%                     somme=somme+2;
%                 end
%             else
%                 if(comp2(2)==comp2(5)) %periodic segment detected
%                     for k3=1:6
%                         if((k3~=2)&&(k3~=5))
%                             if(comp1(k3)==comp2(k3))
%                                 somme=somme+1;
%                             end
%                         end
%                     end
%                     if(comp1(2)==comp1(5))
%                         if(comp1(2)==comp2(2))
%                             somme=somme+2;
%                         end
%                     end
%                     if((comp1(5)==-1)||(comp1(2)==-1))
%                         somme=somme+2;
%                     end
%                 else
                    for k3=1:6
                        if(comp1(k3)==comp2(k3))
                            somme=somme+1;
                        else
                            if((k3==2)||(k3==5))
                                if((comp1(k3)==-1)||(comp2(k3)==-1))
                                    somme=somme+1;
                                    %joker
                                end
                            end
                        end
                    end
%                 end
%             end
            if(somme==6)
                list2=[list2,k2];
            end
        end
    end
    
    sens=[ones(size(list1)),2*ones(size(list2))];
    list1=[list1,list2];
    d=zeros(size(list1));
    for m=1:length(list1)
        d(m)=cotes(k).p1.coord(1)-cotes(list1(m)).p1.coord(1);
    end
    a=find(d<0);
    list1=list1(a);
    sens=sens(a);
    d=abs(d(a));
    b=find(d==min(d))
    list1=list1(b);
    sens=sens(b);
    if(length(list1)>1)
        %choose sens
        sum=0;
        vx=cotes(k).R(end)-cotes(k).R(1);
        vy=cotes(k).Z(end)-cotes(k).Z(1);
        for k5=1:length(cotes(k).R)
            vx1=cotes(k).R(k5)-cotes(k).R(1);
            vy1=cotes(k).Z(k5)-cotes(k).Z(1);
            sum=sum+vx*vy1-vy*vx1;
        end
        sum1=0;
        if(sens(1)==1)
            vx=cotes(list1(1)).R(end)-cotes(list1(1)).R(1);
            vy=cotes(list1(1)).Z(end)-cotes(list1(1)).Z(1);
            for k5=1:length(cotes(list1(1)).R)
                vx1=cotes(list1(1)).R(k5)-cotes(list1(1)).R(1);
                vy1=cotes(list1(1)).Z(k5)-cotes(list1(1)).Z(1);
                sum1=sum1+vx*vy1-vy*vx1;
            end
        else
            vx=cotes(list1(1)).R(1)-cotes(list1(1)).R(end);
            vy=cotes(list1(1)).Z(1)-cotes(list1(1)).Z(end);
            for k5=1:length(cotes(list1(1)).R)
                vx1=cotes(list1(1)).R(k5)-cotes(list1(1)).R(end);
                vy1=cotes(list1(1)).Z(k5)-cotes(list1(1)).Z(end);
                sum1=sum1+vx*vy1-vy*vx1;
            end
        end
        sum2=0;
        if(sens(2)==1)
            vx=cotes(list1(2)).R(end)-cotes(list1(2)).R(1);
            vy=cotes(list1(2)).Z(end)-cotes(list1(2)).Z(1);
            for k5=1:length(cotes(list1(2)).R)
                vx1=cotes(list1(2)).R(k5)-cotes(list1(2)).R(1);
                vy1=cotes(list1(2)).Z(k5)-cotes(list1(2)).Z(1);
                sum2=sum2+vx*vy1-vy*vx1;
            end
        else
            vx=cotes(list1(2)).R(1)-cotes(list1(2)).R(end);
            vy=cotes(list1(2)).Z(1)-cotes(list1(2)).Z(end);
            for k5=1:length(cotes(list1(2)).R)
                vx1=cotes(list1(2)).R(k5)-cotes(list1(2)).R(end);
                vy1=cotes(list1(2)).Z(k5)-cotes(list1(2)).Z(end);
                sum2=sum2+vx*vy1-vy*vx1;
            end
        end
        if(sum~=0)
            if((sum*sum1>0)&&(sum*sum2<0))
                choice=1;
            else
                if((sum*sum1<0)&&(sum*sum2>0))
                    choice=2;
                else
                    %take closest one
                    A=abs([sum-sum1 sum-sum2]);
                    choice=find(A==min(A));
                end
            end
        else
            if(sum1==0)
                choice=1;
            else
                choice=2;
            end
        end
        %         if(sign(sum)==sign(sum1))
        %             list1=list1(1);
        %             sens=sens(1);
        %         else
        %             list1=list1(2);
        %             sens=sens(2);
        %         end
        list1=list1(choice);
        sens=sens(choice);
    end
    for m=1:length(list1);
        
        zones.num=zones.num+1;
        zones.zone(zones.num).coord=[0,0,0,0];
        zones.zone(zones.num).coord(2)=k; %south
        zones.zone(zones.num).coord(1)=list1(m); %north
        zones.zone(zones.num).pA.coord=cotes(k).p1.coord;
        zones.zone(zones.num).pD.coord=cotes(k).p2.coord;
        zones.zone(zones.num).south.R=cotes(k).R;
        zones.zone(zones.num).south.Z=cotes(k).Z;
        if(sens==1)
            zones.zone(zones.num).pB.coord=cotes(list1(m)).p1.coord;
            zones.zone(zones.num).pC.coord=cotes(list1(m)).p2.coord;
            zones.zone(zones.num).north.R=cotes(list1(m)).R;
            zones.zone(zones.num).north.Z=cotes(list1(m)).Z;
        else
            zones.zone(zones.num).pC.coord=cotes(list1(m)).p1.coord;
            zones.zone(zones.num).pB.coord=cotes(list1(m)).p2.coord;
            zones.zone(zones.num).north.R=cotes(list1(m)).R(end:-1:1);
            zones.zone(zones.num).north.Z=cotes(list1(m)).Z(end:-1:1);
        end
        %                 figure(1)
        %                 hold on
        %                 plot(zones.zone(zones.num).south.R,zones.zone(zones.num).south.Z,'r-','LineWidth',2);
        %                 plot(zones.zone(zones.num).north.R,zones.zone(zones.num).north.Z,'b--','LineWidth',2);
        %                 R=(zones.zone(zones.num).south.R(floor(end/2))+zones.zone(zones.num).north.R(floor(end/2)))*0.5;
        %                 Z=(zones.zone(zones.num).south.Z(floor(end/2))+zones.zone(zones.num).north.Z(floor(end/2)))*0.5;
        %                 text(R,Z,num2str(zones.num));
    end
    
    
end

for k=1:zones.num
    %complete
    %west
    comp1=[zones.zone(k).pA.coord,zones.zone(k).pB.coord];
    for n=1:Nseg
        comp2=[segments(n).psimin,segments(n).nX,segments(n).nB,segments(n).type,...
            segments(n).psimax,segments(n).nX,segments(n).nB,segments(n).type];
        somme=0;
        for k3=1:8
            if(comp1(k3)==comp2(k3))
                somme=somme+1;
            else
                if((k3==3)||(k3==7))
                    if((comp1(k3)==-1)||(comp2(k3)==-1))
                        somme=somme+1;
                        %joker
                    end
                end
            end
        end
        if(somme==8)
            zones.zone(k).coord(4)=n;
            zones.zone(k).west.R=segments(n).R;
            zones.zone(k).west.Z=segments(n).Z;
        end
        
    end
    
    
    %complete
    %east
    comp1=[zones.zone(k).pD.coord,zones.zone(k).pC.coord];
    for n=1:Nseg
        comp2=[segments(n).psimin,segments(n).nX,segments(n).nB,segments(n).type,...
            segments(n).psimax,segments(n).nX,segments(n).nB,segments(n).type];
        somme=0;
        
        for k3=1:8
            if(comp1(k3)==comp2(k3))
                somme=somme+1;
            else
                if((k3==3)||(k3==7))
                    if((comp1(k3)==-1)||(comp2(k3)==-1))
                        somme=somme+1;
                        %joker
                    end
                end
            end
        end
        
        if(somme==8)
            zones.zone(k).coord(3)=n;
            zones.zone(k).east.R=segments(n).R;
            zones.zone(k).east.Z=segments(n).Z;
        end
    end
    
    %     plot(zones.zone(k).west.R,zones.zone(k).west.Z,'g-','LineWidth',2);
    %     plot(zones.zone(k).east.R,zones.zone(k).east.Z,'k--','LineWidth',2);
end

%complete cote
for k=1:zones.num
    if(zones.zone(k).pA.coord(4)==1)  %Xpoint
        nX=zones.zone(k).pA.coord(2);
        if(zones.zone(k).pA.coord(1)==X_points.psi(nX))
            zones.zone(k).south.R=[X_points.R(nX),zones.zone(k).south.R];
            zones.zone(k).south.Z=[X_points.Z(nX),zones.zone(k).south.Z];
        end
    end
    if(zones.zone(k).pB.coord(4)==1)  %Xpoint
        nX=zones.zone(k).pB.coord(2);
        if(zones.zone(k).pB.coord(1)==X_points.psi(nX))
            zones.zone(k).north.R=[X_points.R(nX),zones.zone(k).north.R];
            zones.zone(k).north.Z=[X_points.Z(nX),zones.zone(k).north.Z];
        end
    end
    if(zones.zone(k).pC.coord(4)==1)  %Xpoint
        nX=zones.zone(k).pC.coord(2);
        if(zones.zone(k).pC.coord(1)==X_points.psi(nX))
            zones.zone(k).north.R=[zones.zone(k).north.R,X_points.R(nX)];
            zones.zone(k).north.Z=[zones.zone(k).north.Z,X_points.Z(nX)];
        end
    end
    if(zones.zone(k).pD.coord(4)==1)  %Xpoint
        nX=zones.zone(k).pD.coord(2);
        if(zones.zone(k).pD.coord(1)==X_points.psi(nX))
            zones.zone(k).south.R=[zones.zone(k).south.R,X_points.R(nX)];
            zones.zone(k).south.Z=[zones.zone(k).south.Z,X_points.Z(nX)];
        end
    end
end

for k=1:zones.num
    zones.zone(k).south.ismeshed=false;
    zones.zone(k).north.ismeshed=false;
    zones.zone(k).east.ismeshed=false;
    zones.zone(k).west.ismeshed=false;
    zones.zone(k).Xtype_east=segments(zones.zone(k).coord(3)).Xtype;
    zones.zone(k).Xtype_west=segments(zones.zone(k).coord(4)).Xtype;
end

hd=msgbox(['Domain decomposed in ',num2str(zones.num),' zones.']);
WinOnTop(hd);