function [C] = kkb2d(x; y; F; m; n) 
%Berekend de coefficientenmatrix c
%x = [x1...xM] en y = [y1...yN], cartesisch product geeft puntenrooster
% F matrix en m,n veeltermgraden
syms x2 y2
ListXPowers=(x2*ones(1,m+1)).^(0:m);%Generate list of 1, x2,x2^2...till x2^m
ListYPowers=(y2*ones(1,n+1)).^(0:n);%Generate list of 1, y2,y2^2...till y2^n
M=length(x);
N=length(y);
A=zeros(M,m);%Preallocate for speed
B=zeros(N,n);%Preallocate for speed
for i=1:M %For loop adds row per row to matrix A with xi values calculated
    x2=x(i);%Temp value to evaluate list of powers
    A(i,:)=subs(ListXPowers);
end
for i=1:N %For loop adds row per row to matrix B with yi values calculated
    y2=y(i);%Temp value to evaluate list of powers
    B(i,:)=subs(ListYPowers);
end
C=pinv(A)*F*pinv(B); %pinv() gives pseudoinverse of matrix

%Still have to do some shitty tests
end