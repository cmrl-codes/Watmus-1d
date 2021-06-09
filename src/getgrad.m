function [s,dalphadt,dalphan_dalpha0,dalphan_dc]=getgrad(alpha,props,e,de_dc)
[ex,m,adot,h,T,len,nfine,ndf1,dt,dtau,eptol,gtol,bprime]=getprops(props);%return values stored in property array
ndf=size(de_dc,2);
ep(1)=alpha(1);
g(1)=alpha(2);
tau=[0:dtau:T]';
n=size(tau,1);%number of fine scale time points in a cycle
dalphan_dc(1:2,1:ndf,1)=0.0;
dalphan_dalpha0(1:2,1:2,1)=eye(2);
dalph_de=[0 0]';
for i=1:n-1%integration of ep,g and s in a cycle
    epold=ep(i);
    g(i+1)=g(i);
    bkst(i)=0;
    eptr=ep(i);       
    ep(i+1)=finescaleinteg(eptr,ep(i),e(i+1),g(i),bkst(i),props); %backward euler integration of ep 
    dep(i+1)=ep(i+1)-ep(i);
    g(i+1)=g(i)+h*abs(dep(i+1));%integration of g 
    
    sigma=ex*(e(i+1)-ep(i+1));
    
    b=adot*dtau/m*(abs(sigma/g(i+1)))^(1/m-1);
    dalph_dalph(1,1)=1.0+b*ex/g(i+1);
    dalph_dalph(1,2)=b*sigma/g(i+1)^2;
    dalph_dalph(2,1)=-h*sign(sigma);
    dalph_dalph(2,2)=1.0;
    dalph_dalph=inv(dalph_dalph);
    
    temp=eye(2);
    temp(2,1)=-h*sign(sigma);
    
    dalph_dalph_n=dalph_dalph*temp;
    
    dalph_de(1)=b*ex/g(i+1);
    dalph_de(2)=0.0;
    
    dalph_de=dalph_dalph*dalph_de;
    
    dalphan_dc(:,:,i+1)=dalph_de*de_dc(i+1,:)+dalph_dalph_n*dalphan_dc(:,:,i); 
    dalphan_dalpha0(1:2,1:2,i+1)=dalph_dalph_n*dalphan_dalpha0(1:2,1:2,i);
      
end
   
s=ex*(e-ep');%oscillatory stress

dalphadt(1)=ep(n)-ep(1);%change of ep0 in a cycle
dalphadt(2)=g(n)-g(1);%change of g0 in a cycle
end

