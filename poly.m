function [c] = poly(a, T)
%Berekend de klassieke veeltermcoefficienten bij veelterm yn(x).
% Input:
%   a:vector met coefficienten behorende bij stijgende graad van lengte k.
%   T:matrix met coefficienten van Chebyshev-veeltermen van graad nul tot n

%Controleren of de input geldig is
if size(T,1) - length(a) <0
    disp('ERROR: niet genoeg cheb veeltermen voor de a-vector')
    return;
end

%a aanpassen zodat de lengte klopt.
a = [a, zeros(1, size(T,1) - length(a))];

%Preallocate voor snelheid.
MatrixWaarden=zeros(size(a,2));

%Schaal de cheyshev waarden in matrix T met vector a: element tot rij gewijs.
for i =1:size(a,2)
   MatrixWaarden(i,:)=a(i)*T(i,:);
end
%Som alle kolommen om de correcte yn coefficienten te krijgen in stijgende graad.
%Aanname dat sum de som neemt van klein naar groot om de minste fout te krijgen.
c=sum(MatrixWaarden,1);
end