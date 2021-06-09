function coeff=daubdec(f)%Forward transformation using Daubechies-4 wavelets

n=length(f);
coeff(1:n)=0;

nresmin=2;%Two scaling and two wavelets at the coarsest resolution
ftmp=f;

while(log2(n)>=nresmin)
   dbmat_dec=initdb4(n);
   ctmp=dbmat_dec*ftmp(1:n);
   nprev=n;
   n=n/2;
   ftmp(1:n)=ctmp([1:2:nprev]);
   coeff(1:n)=ftmp(1:n);
   coeff(n+1:nprev)=ctmp([2:2:nprev]);
end

end