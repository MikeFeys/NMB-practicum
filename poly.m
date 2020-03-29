%input:
%a-vector from zero to highest degree

function [c] = poly(a, T)
% controle valid input
if size(T,1) - length(a) <0
    disp('ERROR: niet genoeg cheb veeltermen voor de a-vector')
    return;
end
%a is vect with k values
a = [a, zeros(1, size(T,1) - length(a))];

MatrixWaarden=zeros(size(a,2));%Preallocate for speed.
for i =1:size(a,2)
   MatrixWaarden(i,:)=a(i)*T(i,:);%Scale the chebychev values in matrix T with vector a: element to row wise.
end
c=sum(MatrixWaarden,1);%Sum all columns to get the correct yn values. Coef in ascending order of X
end