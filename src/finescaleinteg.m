function y=finescaleinteg(eptr,epold,e,g,bkst,props)
[ex,m,ep0,h,T,len,nfine,ndf,dt,dtau]=getprops(props);
y=fzero(@funcfzero,eptr);
err=funcfzero(y);
    function y=funcfzero(eptr)
        g1=g+h*abs(eptr-epold);
        y=eptr-epold-(ep0*dtau*(((ex*abs((e-eptr))-bkst)/g1)^(1/m))*sign(e-eptr));
    end
end