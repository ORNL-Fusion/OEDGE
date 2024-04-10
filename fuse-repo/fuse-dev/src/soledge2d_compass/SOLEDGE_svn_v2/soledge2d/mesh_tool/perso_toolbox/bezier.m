function bezier=bezier(control)
p=control;
n=length(control);
n1=length(control)-1;
for    i=0:1:n1
sigma(i+1)=factorial(n1)/(factorial(i)*factorial(n1-i));  
end
l=[];
UB=[];
for u=0:0.002:1
for d=1:n
UB(d)=sigma(d)*((1-u)^(n-d))*(u^(d-1));
end
l=cat(1,l,UB); 
end
bezier=l*p;

% figure
% hold on
% plot(bezier(:,1),bezier(:,2),'r')
% plot(p(:,1),p(:,2),'b')
