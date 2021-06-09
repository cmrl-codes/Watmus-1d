function dmat=initdb4mat(nres)%Daubechies-4 forward transformation matrix 

n=2^nres;%Number of data points
dmat(1:n,1:n)=0;
unit(1:n)=0;

for i=1:n
    unit(1:n)=0;
    unit(i)=1;%Unit pulse
    coeff=daubdec(unit');
    dmat(:,i)=coeff;
end
  
end
