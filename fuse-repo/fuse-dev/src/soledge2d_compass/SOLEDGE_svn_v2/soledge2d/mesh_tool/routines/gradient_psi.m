function [x,y]=gradient_psi(pt,psiout,r2D,z2D,flux2D)

    [x1,y1]=follow_grad(pt(1),pt(2),psiout,r2D,z2D,flux2D,3.e-2);
    [x2,y2]=follow_grad(pt(1),pt(2),-10.,r2D,z2D,flux2D,3.e-2);
    x=[x2(end:-1:2),x1];
    y=[y2(end:-1:2),y1];
end