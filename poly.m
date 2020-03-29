function [c] = poly(a, T)
%a is vect with k values
MatrixWaarden=zeros(size(a,2));%Preallocate for speed.
for i =1:size(a,2)
   MatrixWaarden(i,:)=a(i)*T(i,:);%Scale the chebychev values in matrix T with vector a: element to row wise.
end
c=sum(MatrixWaarden,1);%Sum all columns to get the correct yn values.
end