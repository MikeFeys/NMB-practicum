function [e,res,eall]=qr_zonder(A)

% function [e]=qr_zonder(A)
%
% berekent met de QR-methode zonder shift een eigenwaarde van de matrix A
%
% invoer
% A - matrix
% 
% uitvoer
% e - de berekende eigenwaarde
% res - de normen van de residu's voor iedere iteratiestap

[n,m] = size(A);
if n~=m,
  disp('A is geen vierkante matrix')
  return
end
if n<2
  disp('A moet minstens dimensie 2 hebben')
  return
end

res = [];
eall = [];

while any(abs(diag(A,1))>1.e-13)
   res = [res abs(diag(A, 1))];
   [q,r]=qr(A);
   A = r*q;
   eall = [eall, diag(A)];
end
res = [res abs(diag(A,1))];
disp(sprintf('residu = %.1e', abs(A(n,n-1))))
e = diag(A);
[sort(-diag(r)) sort(eig(A))]