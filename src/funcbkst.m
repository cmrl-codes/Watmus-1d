function [y]=func(eptr,epold,e,g,ex,m,ep0,dtau,bkst)%Root solving to obtain plastic strain

y=fzero(@funcfzero,eptr);

    function y=funcfzero(eptr)
    y=eptr-epold-(ep0*dtau*(((ex*abs((e-eptr))-bkst)/g)^(1/m))*sign(e-eptr));
    end
    
end
