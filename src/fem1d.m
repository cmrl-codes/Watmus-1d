clear all;
format long;

nx=3;%Number of nodes
nelx=nx-1;%Number of elements

ex(1:nelx)=200e3;%Elastic modulus of each element
g0(1:nelx)=320;%Yield stress of first element
g0(2)=600;%Yield stress of second element

len(1:nelx)=2.0;%Length of the elements
ar(1:nelx)=1;%Area of the elements

ep0=0.0023;%Plastic multiplier or reference plastic strain rate
h=100;%Hardening rate
xm=0.02;%Exponent

T=2;%Time period

ncyc=1000;%Number of cycles to intergrate
tott=ncyc*T;%Total time - ncyc * T

nres=7;%Number of time points in a cycle (2^nres)
dtau=T/(2^nres-1);%Time step size
tau=[0:dtau:tott]';%Array containing all the time points
nstep=size(tau,1);%Total number of time points

u0=0.008;%Amplitude of applied displacement
omega=2*pi/T;%Angular frequency

u(1:nx,1:nstep)=0;%Initialize nodal displacement array
u(1,:)=0;%Node 1 - BC
u(nx,:)=u0*sin(omega*tau); %Node 3 - BC


ep(1:nelx,1:nstep)=0;%Initialize plastic strain array for all time points
ge(1:nelx,1)=g0(1:nelx);%Initialize yield stress array for first time point

fresid(1:nx)=0;%Initialize residual vector
jacob=zeros(nx,nx);%Initialize Jacobian matrix

maxiter=10; %Maximum number of iterations allowed in a time point
niter=0;%Initialize number of iterations 
err=0;

tic

for n=2:nstep%Time marching
    tt=tau(n);
    u(2:nx-1,n)=u(2:nx-1,n-1);
    iconv=0;
    niter=0;
    while(iconv==0 && niter<maxiter) %NR iteration
       
        fresid(1:nx)=0;
        jacob=zeros(nx,nx);
        for nel=1:nelx
            u1=u(nel,n);
            u2=u(nel+1,n);
            eps=(u2-u1)/len(nel);%Strain in an element 
            
            ep_prev=ep(nel,n-1);%Trial plastic strain
            g=ge(nel,n-1);%Trial yield stress
            ep(nel,n)=funcbkst(ep_prev,ep_prev,eps,g,ex(nel),xm,ep0,dtau,0);
            
            stress=ex(nel)*(eps-ep(nel,n));%Stress
            strel(nel,n)=stress;%Stress in element
            g=g+h*abs(ep(nel,n)-ep_prev);%Update of yield stress 
            ge(nel,n)=g;
            
            dff=stress*ar(nel);
            fresid(nel)=fresid(nel)-dff;%Internal force on local node - 1
            fresid(nel+1)=fresid(nel+1)+dff;%Internal force on local node - 2
            
            b=ep0*dtau/xm*abs(stress/g)^(1/xm-1);
            c=(1+b*h*abs(stress))/g^2;
            dsde=ex/(1+ex*b/(g*(1+c)));%Tangent modulus
            jac=dsde/len(nel);
            
            %Assembly of elemental Jacobian
            jacob(nel,nel)=jacob(nel,nel)+jac;
            jacob(nel,nel+1)=jacob(nel,nel+1)-jac;
            jacob(nel+1,nel)=jacob(nel+1,nel)-jac;
            jacob(nel+1,nel+1)=jacob(nel+1,nel+1)+jac; 
        end
        err=max(abs(fresid(2:nx-1)));%Residual norm
        tol=1e-7;%Tolerance 
        if(err<tol)%Convergence check
            iconv=1;
        else
            niter=niter+1;
            fresid(1)=0;
            fresid(nx)=0;
            jacob(1,1:nx)=0;
            jacob(1:nx,1)=0;
            jacob(nx,1:nx)=0;
            jacob(1:nx,nx)=0;
            jacob(1,1)=1;
            jacob(nx,nx)=1;
            du=inv(jacob)*fresid';%Displacement increment
            u(:,n)=u(:,n)-du(:);%Update of nodal displacement
        end
    end
end

tcpuf = toc;
save single_time_scale_data u ep ge strel tau