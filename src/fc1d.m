clear all;
format long;
load('single_time_scale_data');

global ex;%elastic modulus global declaration
nx=3;%number of nodes
nelx=nx-1;%number of elements
ex=200e3;%elastic modulus
len=2.0;%length of elements
adot=0.0023;%plastic multiplier
h=100;%hardening rate
xm=0.02;%exponent
T=2;%time period
maxiter=10;%maximum number of iterations

u0=0.008;%amplitude of displacement
omega=2*pi/T;%angular frequency

%initializing property array
props(1)=ex;
props(2)=xm;
props(3)=adot;
props(4)=h;
props(5)=T;
props(6)=len;

nres=7;%number of resolution
dt0=1;%cycle jump to commence WATMUS integration
dtau=T/(2^nres-1);%time step size 
tauc=[0:dtau:T]';%time array over a cycle

dt=dt0;
rstep=1;%ratio of current and previous cycle jumps
dmat_full=initdb4mat(nres);%initialize forward transformation matrix


ncycb=2;%Starting cycle
ncyce=1000;%End cycle

%2 previous cycles needed for 2nd order backward difference
t(1)=ncycb;
t(2)=ncycb+dt0;

nfine=2^nres;%Number of time points in a cycle
props(7)=nfine;%Save as properties
nfine1=nfine-1;%Number of increments

ep0st(1:2,1:2)=ep(1:2,nfine1*(t(1:2)-1)+1);%Values of ep0 from previous cycles
g0st(1:2,1:2)=ge(1:2,nfine1*(t(1:2)-1)+1);%Values of g0 from previous cycles

alpha_st(1,1:2,1:2)=ep0st;%Index - 1, ep0, Index - 2, elements, Index - 3, cycles 
alpha_st(2,1:2,1:2)=g0st;%Index - 1, g0, Index - 2, elements, Index - 3, cycles


dep0dt(1:2)=(ep(1:2,nfine1*t(1)+1)-ep(1:2,nfine1*(t(1)-1)+1));%Change of ep0 in a cycle
dg0dt(1:2)=(ge(1:2,nfine1*t(1)+1)-ge(1:2,nfine1*(t(1)-1)+1));%Change of g0 in a cycle

dalphadt_st(1:2,1:2,1)=0;
dalphadt_st(1,1:2,1)=dep0dt(1:2);%Save rate of change of initial value of states
dalphadt_st(2,1:2,1)=dg0dt(1:2);

dep0dt(1:2)=(ep(1:2,nfine1*t(2)+1)-ep(1:2,nfine1*(t(2)-1)+1));
dg0dt(1:2)=(ge(1:2,nfine1*t(2)+1)-ge(1:2,nfine1*(t(2)-1)+1));

dalphadt_st(1,1:2,2)=dep0dt(1:2);
dalphadt_st(2,1:2,2)=dg0dt(1:2);

uc(1,1:nfine)=0.0;%Oscillatory displacement of node 1 for cycle point 2
uc(2,1:nfine)=u(2,nfine1*(t(2)-1)+1:nfine1*t(2)+1);%Oscillatory displacement of node 2 for cycle point 2


u3=u(3,nfine1*(t(2)-1)+1:nfine1*t(2)+1);%Oscillatory displacement of node 3 for cycle point 2
ctemp=dmat_full*uc(2,:)';%Wavelt coefficients of displacement of node 2 for cycle point 2

ctol=1d-4;%Tolerance to select evolving coefficients

ndf=0;
nrem=0;
crem=ctemp;

ftol=1e-7;

props(10)=dtau;%Fine time increment (same as single time scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Selection of evolving coefficients
for i=1:nfine
    if(abs(ctemp(i))>ctol)
        ndf=ndf+1;
        cc(ndf)=ctemp(i);
        dmat(ndf,1:nfine)=dmat_full(i,1:nfine);
        crem(i)=0;
    end
end

c(1,1:ndf)=0.0;
c(2,1:ndf)=cc(1:ndf);
uadd=dmat_full'*crem;%Displacement profile held constant and superposed on evolving profile
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nstep=3;
istop=0;
props(8)=ndf;%number of evolving coefficients
fresid(1:ndf)=0;%initialize residual vector to solve for evolving displacement coefficients of node 2
fresid=fresid';
tic
while(istop==0)
   
   csol=c(nstep-1,1:ndf)';
   niter=0;%initialize number of iterations
   iconv=0;%initialize flag for convergence
   
   while(iconv==0)
       fresid(1:ndf)=0;%initialize
       jac(1:ndf,1:ndf)=0;%initialize
       for nel=1:nelx
          utemp=dmat'*csol+uadd; %displacement profile from evolving and non-evolving coefficients
          alpha_tr(1:2)=alpha_st(:,nel,nstep-1)+dalphadt_st(:,nel,nstep-1)*dt;%trial values of ep0 and g0 for cyclic integration
          eptol=max(1d-8*abs(alpha_st(1,nel,nstep-1)),1d-10);%tolerance for integrating ep0
          gtol=max(1d-6*abs(alpha_st(2,nel,nstep-1)),1d-10);%tolerance integrating g0
          props(11)=eptol;
          props(12)=gtol;
          
          r2=(rstep+1.d0)^2;%factors for backward difference intergration
          r3=r2-(rstep+1);%factors for backward difference intergration

          bprime=(r2-1)/r3;%factors for backward difference intergration
          a1=r2/r3;%factors for backward difference intergration
          a2=1/r3;%factors for backward difference intergration
          
          props(13)=bprime;
          
          %part of backward difference formula having contribution of
          %previous 2 cycle points
          bdf2_prev(1:2)=-a1*alpha_st(:,nel,nstep-1)+a2*alpha_st(:,nel,nstep-2); 
          if(nel==1)
              ef=1.0/len;
              efs=1.0;
              estr=utemp/len;%oscillatory strain in first element
          else
              ef=-1.0/len;
              efs=-1.0;
              estr=(u3'-utemp)/len;%oscillatory strain in second element
          end
          props(9)=dt;%Cycle jump
          
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          %cycle scale integration of ep0 and g0
          %oscillatory stress 
          %derivative of stress with respect to coefficients
          [s,ds_dc,alpha,dalphadt]=constitutive(props,alpha_tr,estr,bdf2_prev,ef*dmat');
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          
          fresid=fresid+efs*dmat*s;%cycle scale residual
          jac=jac+efs*dmat*ds_dc;%cycle scale jacobian
          alpha_st(1:2,nel,nstep)=alpha(1:2);
          dalphadt_st(1:2,nel,nstep)=dalphadt(1:2);
       end
       fmax=max(abs(fresid));%check for convergence
       if(fmax>ftol && niter < maxiter)
           csol=csol-inv(jac)*fresid;
           niter=niter+1;
       else
           iconv=1;
       end
       if niter == maxiter
           iconv=0;
       end
   end
   if(iconv==0)
       inotcon=1
       break
   end
   c(nstep,1:ndf)=csol;
   t(nstep)=t(nstep-1)+dt;
  
   d2epdt2(2,1:2)= (dalphadt_st(1,1:2,nstep)-dalphadt_st(1,1:2,nstep-1))/dt;%second derivative of ep0 wrt cycle at N  
   d2epdt2(1,1:2)= (dalphadt_st(1,1:2,nstep-1)-dalphadt_st(1,1:2,nstep-2))/(t(nstep-1)-t(nstep-2));%second derivative of ep0 wrt cycle at N-DN
   d3epdt3=(d2epdt2(2,1:2)-d2epdt2(1,1:2))/dt;%third derivative of ep0 wrt cycle at N 
   errfac=abs(((rstep+1.d0)^2-(rstep+1.d0)^3)/6.d0);
   d3epdt3=errfac*d3epdt3;
   ep0max=max(max(abs(alpha_st(1,:,nstep))));
   eta=1e-4*ep0max;
   temp=max(abs(d3epdt3(1:2)));%truncation error in backward difference formula
   if(temp>1e-10)
    dtnew=(eta/temp)^(1.0/3.0);%allowable time step (cycle jump)
   else
    dtnew=10000;
   end
   dtnew=floor(dtnew);
   dtmax=min(50,2*dt);%maximum allowed cycle jump
  
   dtmin=dt0;%minimum cycle jump
   if(dtnew>dtmax)
    dtnew=dtmax;
   end
   if(dtnew<dtmin)
    dtnew=dtmin;
   end
   dtprev=dt;
   dt=dtnew;
   rstep=dtprev/dt;  
   if(t(nstep)+dt>ncyce)
       dt=ncyce-t(nstep);
   end
   if(abs(t(nstep)-ncyce)<1e-16)
       istop=1;%stop if end cycle specified is reached
   end
   nstep=nstep+1;%next time point
  
end
tcpuc=toc;

save cscale_data t c alpha_st nres nfine tauc dmat uadd









