clear all

load('single_time_scale_data');

nel=1;%Element number (reinitialize below as well)
var=1;%1 for ep0, 2 for g0 (reinitialize below as well)

T=2;
nres=7;
nfine=2^nres;
nfine1=nfine-1;

ntot=length(u);

ncycst=2;
ncyce=(ntot-1)/nfine1;

cycf=[ncycst:1:ncyce];

j=0;
for i=ncycst:ncyce
    j=j+1;
    if var == 1
    alpha0f(j)=ep(nel,(i-1)*nfine1+1);
    else
    alpha0f(j)=ge(nel,(i-1)*nfine1+1);
    end
end


plot(cycf,alpha0f)
hold on

clear all

load('cscale_data');

nel=1;
var=1;
 
ncyc=length(t);

alpha0c(1:ncyc)=alpha_st(var,nel,1:ncyc);

plot(t,alpha0c,'ro')

hold off




