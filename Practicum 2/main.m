%% Opgave 4
clear
clc
tic
VolSymmetrischeMatrix=load('matrix.txt');
HessenbergVorm=hess(VolSymmetrischeMatrix);
spy(HessenbergVorm,'b',25)
TijdTerInfoInSeconden= toc;
disp(['De tijd die de berekening nam was: ',num2str(TijdTerInfoInSeconden),' seconden.'])
%% Opgave 5
clear
clc
tic
VolSymmetrischeMatrix=load('matrix.txt');
[EigenwaardeQRZonder,ResiduQRZonder]=qr_zonder(VolSymmetrischeMatrix);
[EigenwaardeQRWilkinson,ResiduQRWilkinson]=qr_shiftwilkinsonall(VolSymmetrischeMatrix);
TijdTerInfoInSeconden= toc;
disp(['De tijd die de berekening nam was: ',num2str(TijdTerInfoInSeconden),' seconden.'])
%% Opgave 6 