function [y2] = pccf(x ,y)
[a,b]=size(x);
L=max(a,b);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
if a==1
y1=[conj(y),conj(y)];
   for m=0:L-1
       y2(m+1)=sum(x.*y1(m+1:m+L));
   end

end
if b==1
   y1=[conj(y);conj(y)];
   for m=0:L-1
       y2(m+1)=sum(x.*y1(m+1:m+L));
   end
end
end

