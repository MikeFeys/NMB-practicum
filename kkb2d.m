function [C] = kkb2d(x, y, F, m, n) 
%Berekend de coefficientenmatrix c
%x = [x1...xM] en y = [y1...yN], cartesisch product geeft puntenrooster
% F matrix en m,n veeltermgraden
M=length(x);
N=length(y);
A=zeros(M,m+1);%Preallocate for speed
B=zeros(N,n+1);%Preallocate for speed
if m == n
    if m == n
        for i=1:M %For loop adds row per row to matrix A with xi values calculated
            for j=0:m
                A(i,j+1)=x(i)^j;
                B(i,j+1)=y(i)^j;
            end
        end
    end
else
    for i=1:M %For loop adds row per row to matrix A with xi values calculated
        for j=0:m
            A(i,j+1)=x(i)^j;
        end
    end
    for i=1:N %For loop adds row per row to matrix B with yi values calculated
        for j=0:n
            B(i,j+1)=y(i)^j;
        end
    end
end

C=pinv(A)*F*transpose(pinv(B)); %pinv() gives pseudoinverse of matrix
%Still have to do some shitty tests
end