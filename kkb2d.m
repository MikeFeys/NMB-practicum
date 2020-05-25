function [C] = kkb2d(x, y, F, m, n) 
%Geeft de coefficientenmatrix C terug
% Input:
%   x = [x1...xM] en y = [y1...yN] waarvan cartesisch product puntenrooster geeft
%   F matrix en m,n veeltermgraden

%Waarden toekennen.
M=length(x);
N=length(y);

%Preallocate voor snelheid
A=zeros(M,m+1);
B=zeros(N,n+1);

%Het if statement zorgt voor een snellere berekening indien bepaalde input waarden gelijk
%zijn.
if m == n && M==N
        %Indexeren over de rijen.
        for i=1:M 
            %Indexeren over de elementen in de rij en bijhorende waarden toekennen vanuit de monomiale basis.
            for j=0:m
                A(i,j+1)=x(i)^j;
                B(i,j+1)=y(i)^j;
            end
        end
else
     %Indexeren over de rijen.
    for i=1:M 
        %Indexeren over de elementen in de rij en bijhorende waarden toekennen vanuit de monomiale basis.
        for j=0:m
            A(i,j+1)=x(i)^j;
        end
    end
     %Indexeren over de rijen.
    for i=1:N 
        %Indexeren over de elementen in de rij en bijhorende waarden toekennen vanuit de monomiale basis.
        for j=0:n
            B(i,j+1)=y(i)^j;
        end
    end
end
%Coefficientenmatrix berekenen.
C=pinv(A)*F*transpose(pinv(B));
disp(['Relatieve fout op coefficiëntenmatrix is ',num2str(100*norm(A*C*transpose(B)-F)/norm(C)),'%'])
end