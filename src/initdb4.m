function db4mat_dec = initdb4(n)%Matrix of low and high pass filters at each resolution

  den=4*sqrt(2);
  h(1)=(1+sqrt(3))/den;
  h(2)=(3+sqrt(3))/den;
  h(3)=(3-sqrt(3))/den;
  h(4)=(1-sqrt(3))/den;

  g(1)=(1-sqrt(3))/den;
  g(2)=(sqrt(3)-3)/den;
  g(3)=(3+sqrt(3))/den;
  g(4)=(-1-sqrt(3))/den;

  db4mat_dec(1:n,1:n)=0;
  
  for i=1:2:n-2
     ist=i-1;
     for j=1:4
        db4mat_dec(i,ist+j)=h(j);
        db4mat_dec(i+1,ist+j)=g(j);            
     end
  end

  db4mat_dec(n-1,1)=h(3);
  db4mat_dec(n-1,2)=h(4);
  db4mat_dec(n-1,n-1)=h(1);
  db4mat_dec(n-1,n)=h(2);

  db4mat_dec(n,1)=g(3);
  db4mat_dec(n,2)=g(4);
  db4mat_dec(n,n-1)=g(1);
  db4mat_dec(n,n)=g(2);
   
end
