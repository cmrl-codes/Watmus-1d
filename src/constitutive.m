function [s,ds_dc,alpha,dalphadt,dalphan_dalpha0]=constitutive(props,alpha_tr,estr,bdf2_prev,de_dc)
[ex,m,adot,h,T,len,nfine,ndf,dt,dtau,eptol,gtol,beta]=getprops(props);%return values stored in property array
ndf=size(de_dc,2);
iconv=0;%initialize convergence flag

gamma=beta+dt;
alpha=alpha_tr;
while(iconv==0)%check for convergence
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %change of ep0 and g0 in a cycle - dy0/dN
    %oscillatory stress
    %derivative of oscillatory ep and g wrt ep0 and g0
    %derivative of oscillatory ep and g wrt to coefficients
    [s,dalphadt,dalphan_dalpha0,dalphan_dc]=getgrad(alpha,props,estr,de_dc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fc=beta*alpha+bdf2_prev-dalphadt*dt;%residual to solve for ep0 and g0
    
    if(abs(fc(1))>=eptol||abs(fc(2))>=gtol) %check with tolerance
        dfc_dalpha0=gamma*eye(2)-dalphan_dalpha0(:,:,nfine)*dt;
        dalpha=-inv(dfc_dalpha0)*fc';
        alpha=alpha+dalpha';%iterative update of ep0 and g0
    else
        iconv=1;
    end
end
if(iconv==0)
    inotcon=1
    pause;
end

fac=dt/gamma;
dalpha0_dc=fac*inv(eye(2)-fac*dalphan_dalpha0(:,:,nfine))*dalphan_dc(:,:,nfine);%derivate of ep0 and g0 wrt coefficients
for i=1:nfine
    temp=dalphan_dalpha0(:,:,i)*dalpha0_dc+dalphan_dc(:,:,i); 
    ds_dc(i,1:ndf)=ex*(de_dc(i,1:ndf)-temp(1,1:ndf));%derivate of stress wrt coefficients     
end
end